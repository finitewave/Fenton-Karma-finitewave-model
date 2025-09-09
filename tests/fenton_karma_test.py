import numpy as np
import pytest

from implementation.fenton_karma_0d import FentonKarma0D, Stimulation


def prepare_model(model_class, dt, curr_dur, curr_value, t_prebeats):
    """
    Prepares a 2D cardiac model with a stimulation protocol.

    Parameters
    ----------
    model_class : Callable
        The cardiac model class to be instantiated.
    dt : float
        Time step for the simulation (ms or model units).
    curr_value : float
        Amplitude of the stimulus current (μA/cm² or model units).
    curr_dur : float
        Duration of each stimulus pulse (ms or model units).
    t_prebeats : float
        Interval between preconditioning stimuli (ms or model units).

    Returns
    -------
    model : Model
        Configured and initialized model ready for simulation.
    """

    stimulations = [Stimulation(t_start=0.0, duration=curr_dur, amplitude=curr_value),
                    Stimulation(t_start=t_prebeats, duration=curr_dur, amplitude=curr_value),
                    Stimulation(t_start=2*t_prebeats, duration=curr_dur, amplitude=curr_value),
                    Stimulation(t_start=3*t_prebeats, duration=curr_dur, amplitude=curr_value)]

    model = model_class(dt=dt, stimulations=stimulations)

    return model


def calculate_apd(u, dt, threshold, beat_index=3):
    """
    Calculates the action potential duration (APD) for a single beat (third by default).

    Parameters
    ----------
    u : np.ndarray
        Membrane potential time series.
    dt : float
        Time step of the simulation (ms).
    threshold : float
        Voltage threshold to define APD90 (e.g., -70 mV or 0.1 for normalized models).
    beat_index : int, optional
        Index of the beat to analyze (default is 3).

    Returns
    -------
    apd : float or None
        Duration of the action potential (ms or model units), or None if no complete AP was found.
    """
    up_idx = np.where((u[:-1] < threshold) & (u[1:] >= threshold))[0]
    down_idx = np.where((u[:-1] > threshold) & (u[1:] <= threshold))[0]

    if len(up_idx) <= beat_index or len(down_idx) == 0:
        return None

    ap_start = up_idx[beat_index]
    ap_end_candidates = down_idx[down_idx > ap_start]
    if len(ap_end_candidates) == 0:
        return None

    ap_end = ap_end_candidates[0]
    return (ap_end - ap_start) * dt


def test_model_attributes():
    """
    Test that the model has the expected attributes.
    Checks for the presence of key variables and parameters in the 0D Model.
    """
    model = FentonKarma0D(dt=0.01, stimulations=[])

    assert 'u' in model.variables, "Fenton-Karma should have variable 'u'"
    assert 'v' in model.variables, "Fenton-Karma should have variable 'v'"
    assert 'w' in model.variables, "Fenton-Karma should have variable 'w'"

    expected_params = ['tau_r', 'tau_o', 'tau_d', 'tau_si', 'tau_v_m',
                          'tau_v_p', 'tau_w_m', 'tau_w_p', 'k', 'u_c', 'uc_si']
     
    for param in expected_params:
        assert param in model.parameters, f"Fenton-Karma should have parameter '{param}'"


def test_model_run():
    """
    Test the model run for a short simulation.
    Runs the 0D Model with a predefined stimulation protocol and checks
    that the membrane potential 'u' stays within expected physiological ranges.
    """
    t_prebeats = 1000.0 # interval between preconditioning stimuli (ms or model units).
    t_calc = 1000.0     # time after the last preconditioning beat to continue recording (ms or model units).
    t_max = 3*t_prebeats + t_calc
    model = prepare_model(FentonKarma0D, dt=0.01, curr_dur=0.2, curr_value=1.0, t_prebeats=t_prebeats)
    model.run(t_max=t_max)
    u = np.array(model.history['u'])

    assert np.max(u) == pytest.approx(1.0, abs=0.1)
    assert np.min(u) == pytest.approx(0.0, abs=0.01)

    apd = calculate_apd(u, model.dt, threshold=0.1)
    assert 200 <= apd <= 300, f"Fenton-Karma is out of expected range {apd}"
