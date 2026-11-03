"""Checks against Williams 2013 Table 1, Eqs. 1-14 and Fig. 5F."""
import math
import numpy as np
import pytest
from fenton_karma import ops
from implementation.fenton_karma_0d import FentonKarma0D, Stimulation


def baseline_rates(V, light, p_chr2):
    pars = ops.get_parameters()
    flux = ops.calc_flux(pars['sigma_ret'], light, pars['lambda_nm'], pars['w_loss'], pars['hc'])
    return (ops.calc_k1(pars['eps1'], flux, p_chr2),
            ops.calc_k2(pars['eps2'], flux, p_chr2), ops.calc_Gd1(V),
            ops.calc_Gd2(), ops.calc_Gr(V), ops.calc_e12(light), ops.calc_e21(light))


def current(g, o1, gamma, o2, V):
    return ops.calc_J_chr2(g, o1, gamma, o2, ops.calc_G_times_V(V))


def test_published_rates_and_temperature_scaling():
    rates = baseline_rates(-80.0, 1.0, 1.0)
    # Independently evaluated Table 1 values, including physical photon units.
    assert rates == pytest.approx((
        0.3147133317296513, 0.05162257345301838, 0.1177873544085294,
        0.05, 0.00023606926407714654, 0.029767089876257535,
        0.02301367190100603,
    ))
    # Infer each transition from pure states to check scaling in ionic_step.
    for temperature in [22.0, 32.0, 37.0]:
        pars = ops.get_parameters()
        pars.update(temperature=temperature, irradiance=1.0)
        x = (temperature - 22.0) / 10.0
        k1, k2, gd1, gd2, gr, e12, e21 = rates
        k1 *= 1.46**x
        k2 *= 2.77**x
        gd1 *= 1.97**x
        gd2 *= 1.77**x
        gr *= 2.56**x
        e12 += .011*(1.1**x-1)
        e21 += .008*(1.95**x-1)
        transitions = [(-k1,k1,0,0), (gd1,-gd1-e12,e12,0),
                       (0,e21,-gd2-e21,gd2), (gr,0,k2,-k2-gr)]
        for occupied, expected in zip(['c1','o1','o2','c2'], transitions):
            state = ops.get_variables()
            state.update(u=.05, p_chr2=1, c1=0, o1=0, o2=0, c2=0)
            state[occupied] = 1
            result = ops.ionic_step(.01, **state, **pars)
            derivative = [(value-state[key])/.01 for key,value in zip(['c1','o1','o2','c2'],result[3:7])]
            assert derivative == pytest.approx(expected, abs=1e-12)
    assert ops.calc_dp_chr2(0.0, 0.1, 1.3) == pytest.approx(0.5 / 1.3)


def test_current_voltage_units_and_zero_voltage_limit():
    for V in [-85.0, -20.0, 20.0]:
        G = (10.6408 - 14.6408 * math.exp(-V / 42.7671)) / V
        assert current(0.4, 0.2, 0.1, 0.3, V) == pytest.approx(0.4 * G * 0.23 * V)
    assert current(0.4, 0.2, 0.1, 0.3, 0.0) == pytest.approx(-0.368)
    state = ops.get_variables()
    state.update(o1=0.2, c1=0.8)
    pars = ops.get_parameters()
    lit = ops.ionic_step(0.01, **state, **pars)
    pars['g_chr2'] = 0.0
    disabled = ops.ionic_step(0.01, **state, **pars)
    assert lit[0] - disabled[0] == pytest.approx(-current(0.4, .2, .1, 0, -85)/100)
    assert lit[1:] == disabled[1:]


def test_dark_channels_preserve_original_fk_rhs_and_gates():
    pars = ops.get_parameters()
    for u in [0.0, 0.05, 0.25, 0.9, 1.0]:
        state = ops.get_variables()
        state.update(u=u, v=0.7, w=0.6)
        result = ops.ionic_step(.01, **state, **pars)
        p, q = float(u >= .25), float(u <= .0025)
        Jfi = -.7*p*(u-.1)*(.97-u)/.175
        Jso = u*(1-p)*.3/4.5+p/(5+15*(1-math.tanh(50*(u-.85))))
        Jsi = 0 if u < .05 else -.6*(1-math.tanh(10*(u-.9)))*(1+math.tanh(7*(u-.35)))/(62*(1+math.exp(4.5*(u-.9))))
        assert result[:3] == pytest.approx((
            -Jfi-Jso-Jsi,
            .7+.01*((1-p)*.3/(40*q+(1-q)*2000)-p*.7/10),
            .6+.01*((1-p)*(1-.6**4)/305-p*.6/320),
        ))
        assert result[3:7] == (1., 0., 0., 0.)


def test_all_transitions_conserve_probability():
    c1, o1, o2, c2 = .4, .2, .1, .3
    k1, k2, gd1, gd2, gr, e12, e21 = baseline_rates(-40, 5.5, .8)
    derivatives = (
        ops.calc_dc1(c1, c2, o1, gr, gd1, k1),
        ops.calc_do1(c1, o1, o2, k1, gd1, e12, e21),
        ops.calc_do2(c2, o1, o2, k2, gd2, e12, e21),
        ops.calc_dc2(c2, o2, k2, gd2, gr),
    )
    assert sum(derivatives) == pytest.approx(0, abs=1e-15)
    # C2 must open through k2, even when neither open state is occupied.
    assert ops.calc_do2(1, 0, 0, k2, gd2, e12, e21) == k2


@pytest.mark.parametrize('temperature', [22., 37.])
def test_optical_pulse_excites_and_channels_close(temperature):
    model = FentonKarma0D(.01, [], [Stimulation(10, 20, 1)])
    model.parameters['temperature'] = temperature
    model.run(600)
    history = model.history
    assert max(history['u']) > .8
    assert history['u'][-1] < .01
    probabilities = np.array([history[k] for k in ['c1', 'o1', 'o2', 'c2']])
    assert np.all(np.isfinite(probabilities))
    assert probabilities.min() >= -1e-12
    assert probabilities.max() <= 1 + 1e-12
    np.testing.assert_allclose(probabilities.sum(axis=0), 1, atol=1e-12)
    assert history['o1'][-1] + history['o2'][-1] < 1e-8


def test_pulse_time_step_convergence():
    peaks = []
    for dt in [.02, .01]:
        model = FentonKarma0D(dt, [], [Stimulation(10, 20, 1)])
        model.run(100)
        peaks.append(max(model.history['u']))
    assert abs(peaks[0] - peaks[1]) < .001


def test_photon_units_and_exported_parameter_contract():
    import inspect
    pars = ops.get_parameters()
    assert ops.calc_flux(pars['sigma_ret'], 1., pars['lambda_nm'], pars['w_loss'], pars['hc']) == pytest.approx(.3687326675215598)
    names = set(inspect.signature(ops.ionic_step).parameters)
    assert names == {'dt'} | set(pars) | set(ops.get_variables())
    for name in ops.__all__:
        assert callable(getattr(ops, name))
    # Changing a physical parameter reaches the state update.
    state = ops.get_variables()
    state['p_chr2'] = 1.
    pars['irradiance'] = 1.
    first = ops.ionic_step(.01, **state, **pars)
    pars['sigma_ret'] *= 2
    second = ops.ionic_step(.01, **state, **pars)
    assert second[4] == pytest.approx(2 * first[4])
