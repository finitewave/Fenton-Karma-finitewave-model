"""
ops.py — mathematical core of the __template__ model.

This module provides functions to compute the __template__ model equations,
as well as functions to retrieve default parameters and initial
values for the state variables.

The __template__ model is a simplified representation of cardiac electrical
activity, capturing essential features of excitation and recovery in cardiac tissue.

References:
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",  # add other calc_* functions as needed
    # e.g. "calc_dv",
)


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.
    Example (ionic model): {"u": -84.0, "m": 0.01, "h": 0.99}
    Example (phenomenological): {"u": 0.0, "v": 0.1}
    """
    raise NotImplementedError("The get_variables method must be implemented in a subclass.")


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.
    Example (ionic model): {"g_Na": 120.0, "E_Na": 50.0, "C_m": 1.0}
    Example (phenomenological): {"a": 0.1, "b": 0.5}
    """
    raise NotImplementedError("The get_parameters method must be implemented in a subclass.")


def calc_rhs() -> float:
    """
    Computes the right-hand side of the __template__ model.
    """
    raise NotImplementedError("The calc_rhs method must be implemented in a subclass.")

