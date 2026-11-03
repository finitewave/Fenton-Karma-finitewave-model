"""Porcine ventricular minimal model with Williams et al. ChR2(H134R).

Voltage u is dimensionless; time is in ms. ChR2 uses V = V_rest + V_scale*u
in mV and contributes -I_chr2/(C_m*V_scale) to du/dt. The voltage mapping
is a configurable coupling assumption, not a fit supplied by either paper.

References:
Kaboudian et al. (2019), doi:10.1126/sciadv.aav6019, pp. 3-4.
Williams et al. (2013), doi:10.1371/journal.pcbi.1003220,
Eqs. 1-14, Table 1 and Fig. 5F. ChR2 parameters default to 22 degrees C.
"""

__all__ = (
    "get_diffusion_coefficient",
    "get_variables",
    "get_parameters",
    "ionic_step",
    "calc_rhs",
    "calc_where",
    "calc_Jfi",
    "calc_Jso",
    "calc_Jsi",
    "calc_dv",
    "calc_dw",
    "calc_tau_v_m",
    "calc_tau_si",
    "calc_J_chr2",
    "calc_G_times_V",
    "calc_Gd1",
    "calc_Gd2",
    "calc_Gr",
    "calc_e12",
    "calc_e21",
    "calc_k1",
    "calc_k2",
    "calc_flux",
    "calc_dc1",
    "calc_do1",
    "calc_do2",
    "calc_dc2",
    "calc_dp_chr2",
)

from math import tanh, exp, log

def get_diffusion_coefficient() -> dict[str, float]:
    return {"D_model": 0.1}

def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.
    """
    return {
        "u": 0.0, "v": 1.0, "w": 1.0,
        "c1": 1.0, "o1": 0.0, "o2": 0.0, "c2": 0.0, "p_chr2": 0.0,
    }


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.
    """
    return {
        "k": 10.0,
        "g_fi": 1 / 0.175,
        "tau_v_p": 10.0,
        "tau_w_p": 320.0,
        "tau_w_m": 305.0,
        "u_si": 0.9,
        "u_p": 0.25,
        "u_q": 0.0025,
        "irradiance": 0.0,  # mW/mm^2; may be varied between time steps
        "g_chr2": 0.4,  # mS/cm^2
        "gamma_chr2": 0.1,
        "sigma_ret": 12e-20,  # m^2, Table 1 (12, not 1.2)
        "lambda_nm": 470.0,  # nm
        "w_loss": 0.77,  # division convention for photon losses
        "hc": 1.986446e-25,  # J*m
        "eps1": 0.8535,  # quantum efficiency at 22 C
        "eps2": 0.14,
        "tau_chr2": 1.3,  # ms
        "temperature": 22.0,  # degrees C; only ChR2 rates are scaled
        "V_rest": -85.0,  # mV, coupling assumption
        "V_scale": 100.0,  # mV per unit u, must be positive
        "C_m": 1.0,  # uF/cm^2, must be positive

    }


def ionic_step(
    dt,
    u,
    v,
    w,
    c1,
    o1,
    o2,
    c2,
    p_chr2,
    k,
    g_fi,
    tau_v_p,
    tau_w_p,
    tau_w_m,
    u_si,
    u_p,
    u_q,
    irradiance,
    g_chr2,
    gamma_chr2,
    temperature,
    sigma_ret,
    lambda_nm,
    w_loss,
    hc,
    eps1,
    eps2,
    tau_chr2,
    V_rest,
    V_scale,
    C_m,
):
    """Return (rhs, v_new, w_new, c1_new, o1_new, o2_new, c2_new,
    p_chr2_new) using forward Euler with all currents evaluated at the old state.

    ChR2 states are probabilities, initially entirely in C1. Irradiance is
    nonnegative in mW/mm^2. V_rest and V_scale map u to mV; C_m is in uF/cm^2.
    The fitted ChR2 reversal parameter is fixed at 0 mV (see calc_J_chr2).
    """
    p = calc_where(u >= u_p, 1.0, 0.0)
    q = calc_where(u <= u_q, 1.0, 0.0)

    tau_v_m = calc_tau_v_m(q)
    tau_si = calc_tau_si(u, u_si)

    J_fi = calc_Jfi(u, v, p, g_fi)
    J_so = calc_Jso(u, v, p)
    J_si = calc_Jsi(u, w, k, u_si, tau_si)

    dv = calc_dv(v, p, tau_v_m, tau_v_p)
    dw = calc_dw(w, p, tau_w_m, tau_w_p)

    V = V_rest + V_scale * u
    irradiance = calc_where(irradiance > 0.0, irradiance, 0.0)
    theta = 100.0 * irradiance  # Table 1, optical activation protocol
    exponent = (temperature - 22.0) / 10.0
    flux = calc_flux(sigma_ret, irradiance, lambda_nm, w_loss, hc)
    k1 = calc_k1(eps1, flux, p_chr2) * 1.46**exponent
    k2 = calc_k2(eps2, flux, p_chr2) * 2.77**exponent
    Gd1 = calc_Gd1(V) * 1.97**exponent
    Gd2 = calc_Gd2() * 1.77**exponent
    Gr = calc_Gr(V) * 2.56**exponent
    e12 = calc_e12(irradiance) + 0.011 * (1.1**exponent - 1.0)
    e21 = calc_e21(irradiance) + 0.008 * (1.95**exponent - 1.0)

    v_new = v + dv * dt
    w_new = w + dw * dt

    c1_new = c1 + dt * calc_dc1(c1, c2, o1, Gr, Gd1, k1)
    o1_new = o1 + dt * calc_do1(c1, o1, o2, k1, Gd1, e12, e21)
    o2_new = o2 + dt * calc_do2(c2, o1, o2, k2, Gd2, e12, e21)
    c2_new = c2 + dt * calc_dc2(c2, o2, k2, Gd2, Gr)
    p_chr2_new = p_chr2 + dt * calc_dp_chr2(p_chr2, theta, tau_chr2)
    G_times_V = calc_G_times_V(V)
    J_chr2 = calc_J_chr2(g_chr2, o1, gamma_chr2, o2, G_times_V)
    rhs = calc_rhs(J_fi, J_so, J_si) - J_chr2 / (C_m * V_scale)

    return rhs, v_new, w_new, c1_new, o1_new, o2_new, c2_new, p_chr2_new


def calc_rhs(J_fi, J_so, J_si) -> float:
    """
    Computes the right-hand side of the Fenton-Karma model.

    Parameters
    ----------
    J_fi : float
        Fast inward current.
    J_so : float
        Slow outward current.
    J_si : float
        Slow inward current.

    Returns
    -------
    float
        The net change in membrane potential.
    """
    return -J_fi - J_so - J_si


def calc_where(cond, x, y):
    if cond:
        return x
    return y


def calc_Jfi(u, v, p, g_fi):
    """
    Computes the fast inward current (J_fi) for the Fenton-Karma model.

    This current is responsible for the rapid depolarization of the membrane
    potential. It is active only when the membrane potential exceeds a threshold `u_c`.

    Parameters
    ----------
    u : float
        Current membrane potential (dimensionless).
    v : float
        Fast recovery gate (sodium channel inactivation).
    p : float
        Indicator of whether the membrane potential is above threshold (1.0) or below (0.0).
    g_fi : float
        Conductance parameter of the fast inward current.

    Returns
    -------
    float
        Value of the fast inward current at this point.
    """
    return -(v * p * (u - 0.1) * (0.97 - u)) * g_fi


def calc_Jso(u, v, p):
    """
    Computes the slow outward current (J_so) for repolarization.

    This current contains two parts:
    - a linear repolarizing component active below threshold `u_c`
    - a constant repolarizing component above threshold

    Parameters
    ----------
    u : float
        Membrane potential.
    v : float
        Slow recovery gate.
    p : float
        Indicator of whether the membrane potential is above threshold (1.0) or below (0.0).

    Returns
    -------
    float
        Value of the outward repolarizing current.
    """

    return u * (1.0 - p) * (1.0 - v) / 4.5 + p / (5 + 15 * (1.0 - tanh(50 * (u - 0.85))))


def calc_Jsi(u, w, k, u_si, tau_si):
    """
    Computes the slow inward (calcium-like) current (J_si).

    This current is responsible for the plateau phase of the action potential
    and depends on the gating variable `w` and a smoothed activation threshold.

    Parameters
    ----------
    u : float
        Membrane potential.
    w : float
        Slow recovery gate.
    k : float
        Steepness of the tanh activation curve.
    u_si : float
        Activation threshold for the slow inward current.
    tau_si : float
        Time constant for the slow inward current.

    Returns
    -------
    float
        Value of the slow inward current.
    """
    J_si = (
        - w
        * (1. - tanh(k * (u - u_si)))
        * (1. + tanh(7 * (u - 0.35)))
        / tau_si
    )
    J_si = calc_where(u < 0.05, 0.0, J_si)
    return J_si


def calc_dv(v, p, tau_v_m, tau_v_p):
    """
    Calculates the fast recovery gate `v`.

    This gate controls sodium channel availability and changes depending on
    whether the membrane potential is below or above a critical threshold.

    Parameters
    ----------
    v : float
        Current value of the recovery variable.
    p : float
        Indicator of whether the membrane potential is above threshold (1.0) or below (0.0).
    tau_v_m : float
        Time constant below threshold.
    tau_v_p : float
        Time constant above threshold.

    Returns
    -------
    float
        Updated value of `v`.
    """
    return (1.0 - p) * (1.0 - v) / tau_v_m - p * v / tau_v_p


def calc_tau_v_m(q):
    return 40. * q + (1 - q) * 2000.


def calc_tau_si(u, u_si):
    return 62 * (1. + exp(4.5 * (u - u_si)))


def calc_dw(w, p, tau_w_m, tau_w_p):
    """
    Calculates the slow recovery gate `w`.

    This gate represents the calcium channel recovery and decays similarly to `v`,
    depending on whether the membrane potential is above or below threshold `u_c`.

    Parameters
    ----------
    w : float
        Current value of the recovery variable.
    p : float
        Indicator of whether the membrane potential is above threshold (1.0) or below (0.0).
    tau_w_m : float
        Time constant below threshold.
    tau_w_p : float
        Time constant above threshold.

    Returns
    -------
    float
        Updated value of `w`.
    """
    return (1.0 - p) * (1.0 - w**4) / tau_w_m - p * w / tau_w_p


def calc_J_chr2(g_chr2, o1, gamma, o2, G_times_V):
    """Eq. 1 current density (uA/cm^2) for the paper's fixed E_chr2=0 mV."""
    return g_chr2 * (o1 + gamma * o2) * G_times_V


def calc_G_times_V(V):
    """Eq. 12 times V, in mV; this is NOT G(V) alone.

    With E_chr2=0, cancel V analytically in G(V)*(V-E_chr2), including
    its finite limit at V=0. The empirical fit crosses zero near +13.65 mV.
    """
    return 10.6408 - 14.6408 * exp(-V / 42.7671)


def calc_Gd1(V):
    return (0.075 + 0.043 * tanh((V + 20.0) / -20.0))


def calc_Gd2():
    return 0.05


def calc_Gr(V):
    return 4.34587e-5 * exp(-0.0211539274 * V)


def calc_e12(light):
    return 0.011 + 0.005 * log(1.0 + light / 0.024)


def calc_e21(light):
    return 0.008 + 0.004 * log(1.0 + light / 0.024)


def calc_k1(eps1, flux, p_chr2):
    return eps1 * flux * p_chr2


def calc_k2(eps2, flux, p_chr2):
    return eps2 * flux * p_chr2


def calc_flux(sigma_ret, light, lambda_nm, w_loss, hc):
    """Eq. 9 photon absorption rate per molecule, in ms^-1.

    Irradiance is mW/mm^2 and wavelength is nm. Irradiance conversion to
    W/m^2 (1e3) cancels conversion from s^-1 to ms^-1 (1e-3);
    wavelength still needs the factor 1e-9 to convert nm to m.
    """
    return sigma_ret * light * lambda_nm * 1e-9 / (w_loss * hc)


def calc_dc1(c1, c2, o1, Gr, Gd1, k1):
    return Gr * c2 + Gd1 * o1 - k1 * c1


def calc_do1(c1, o1, o2, k1, Gd1, e12, e21):
    return k1 * c1 - (Gd1 + e12) * o1 + e21 * o2


def calc_do2(c2, o1, o2, k2, Gd2, e12, e21):
    return k2 * c2 - (Gd2 + e21) * o2 + e12 * o1


def calc_dc2(c2, o2, k2, Gd2, Gr):
    return Gd2 * o2 - (k2 + Gr) * c2


def calc_dp_chr2(p_chr2, theta, tau_chr2):
    """Eqs. 10-11; theta=100*irradiance is supplied by ionic_step."""
    S0 = 0.5 * (1.0 + tanh(120.0 * (theta - 0.1)))
    return (S0 - p_chr2) / tau_chr2
