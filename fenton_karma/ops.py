"""
ops.py — mathematical core of the Fenton-Karma model.

This module provides functions to compute the Fenton-Karma model equations,
as well as functions to retrieve default parameters and initial
values for the state variables.

The Fenton-Karma model is a minimal three-variable model designed to reproduce
essential features of human ventricular action potentials, including restitution, 
conduction velocity dynamics, and spiral wave behavior. It captures the interaction 
between fast depolarization, slow repolarization, and calcium-mediated effects 
through simplified phenomenological equations.

This implementation corresponds to the MLR-I parameter set described in the original paper
and supports 2D isotropic and anisotropic tissue simulations with diffusion.

References:
- Fenton, F. H., & Karma, A. (1998). Vortex dynamics in three-dimensional
    continuous myocardium with fiber rotation: Filament instability and fibrillation.
    Chaos: An Interdisciplinary Journal of Nonlinear Science, 8(1), 20-47.

DOI: https://doi.org/10.1063/1.166311
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",
    "calc_Jfi",
    "calc_Jso",
    "calc_Jsi",
    "calc_dv",
    "calc_dw",
)

from math import tanh


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.
    """
    return {"u": 0.0, "v": 1.0, "w": 1.0}


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.
    """
    return {"tau_r": 130.0, "tau_o": 12.5, "tau_d": 0.172, "tau_si": 127.0,
            "tau_v_m": 18.2, "tau_v_p": 10.0, "tau_w_m": 80.0, "tau_w_p": 1020.0,
            "k": 10.0, "u_c": 0.13, "uc_si": 0.85}


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


def calc_Jfi(u, v, u_c, tau_d):
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
    u_c : float
        Activation threshold for the inward current.
    tau_d : float
        Time constant for depolarization.

    Returns
    -------
    float
        Value of the fast inward current at this point.
    """
    H = 1.0 if (u - u_c) >= 0 else 0.0
    return -(v*H*(1-u)*(u - u_c))/tau_d


def calc_Jso(u, u_c, tau_o, tau_r):
    """
    Computes the slow outward current (J_so) for repolarization.

    This current contains two parts:
    - a linear repolarizing component active below threshold `u_c`
    - a constant repolarizing component above threshold

    Parameters
    ----------
    u : float
        Membrane potential.
    u_c : float
        Activation threshold.
    tau_o : float
        Time constant for subthreshold repolarization.
    tau_r : float
        Time constant for suprathreshold repolarization.

    Returns
    -------
    float
        Value of the outward repolarizing current.
    """
    H1 = 1.0 if (u_c - u) >= 0 else 0.0
    H2 = 1.0 if (u - u_c) >= 0 else 0.0

    return u*H1/tau_o + H2/tau_r


def calc_Jsi(u, w, k, uc_si, tau_si):
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
    uc_si : float
        Activation threshold for the slow inward current.
    tau_si : float
        Time constant for the slow inward current.

    Returns
    -------
    float
        Value of the slow inward current.
    """
    return -w*(1 + tanh(k*(u - uc_si)))/(2*tau_si)


def calc_dv(v, u, u_c, tau_v_m, tau_v_p):
    """
    Calculates the fast recovery gate `v`.

    This gate controls sodium channel availability and changes depending on
    whether the membrane potential is below or above a critical threshold.

    Parameters
    ----------
    v : float
        Current value of the recovery variable.
    u : float
        Membrane potential.
    u_c : float
        Activation threshold.
    tau_v_m : float
        Time constant below threshold.
    tau_v_p : float
        Time constant above threshold.

    Returns
    -------
    float
        Updated value of `v`.
    """
    H1 = 1.0 if (u_c - u) >= 0 else 0.0
    H2 = 1.0 if (u - u_c) >= 0 else 0.0
    return H1*(1 - v)/tau_v_m - H2*v/tau_v_p


def calc_dw(w, u, u_c, tau_w_m, tau_w_p):
    """
    Calculates the slow recovery gate `w`.

    This gate represents the calcium channel recovery and decays similarly to `v`,
    depending on whether the membrane potential is above or below threshold `u_c`.

    Parameters
    ----------
    w : float
        Current value of the recovery variable.
    u : float
        Membrane potential.
    u_c : float
        Activation threshold.
    tau_w_m : float
        Time constant below threshold.
    tau_w_p : float
        Time constant above threshold.

    Returns
    -------
    float
        Updated value of `w`.
    """
    H1 = 1.0 if (u_c - u) >= 0 else 0.0
    H2 = 1.0 if (u - u_c) >= 0 else 0.0
    return H1*(1 - w)/tau_w_m - H2*w/tau_w_p
