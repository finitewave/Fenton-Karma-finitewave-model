"""
This module provides a simple interface to run the Fenton-Karma model in a 0D setting,
i.e., without spatial dimensions. It includes class for defining stimulation protocols
and a class for the 0D model itself.

"""

from fenton_karma import ops


class Stimulation:
    """
    Stimulus protocol for the 0D model.

    Parameters
    ----------
    t_start : float
        Start time (ms) of the first stimulus window.
    duration : float
        Duration (ms) of a single pulse.
    amplitude : float
        Pulse amplitude in the same units as du/dt contribution (typically "units/ms").

    Method
    ------
    stim(t: float) -> float
        Returns the instantaneous stimulus value at time t.

    """

    def __init__(self, t_start: float, duration: float, amplitude: float):
        self.t_start = t_start
        self.duration = duration
        self.amplitude = amplitude

    def stim(self, t: float) -> float:
        return self.amplitude if self.t_start <= t < self.t_start + self.duration else 0.0


class FentonKarma0D:
    """
    Fenton-Karma model in 0D.

    Parameters
    ----------

    dt : float
        Time step size (ms).
    stimulations : list[Stimulation]
        List of stimulation protocols to apply during the simulation.

    Attributes
    ----------
    variables : dict[str, float]
        Current state variables of the model.
    parameters : dict[str, float]
        Model parameters.
    history : dict[str, list[float]]
        Time history of state variables for post-processing.
    
    Methods
    -------
    step(i: int)
        Perform a single time step update.
    run(t_max: float)
        Run the simulation up to time t_max.
    """
    def __init__(self, dt: float, stimulations: list[Stimulation]):
        self.dt = dt
        self.stimulations = stimulations
        self.variables = ops.get_variables()
        self.parameters = ops.get_parameters()
        self.history = {s: [] for s in self.variables}

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """
        self.variables["v"] += self.dt *ops.calc_dv(self.variables["v"], self.variables["u"], 
                                          self.parameters["u_c"], self.parameters["tau_v_m"], self.parameters["tau_v_p"])
        self.variables["w"] += self.dt *ops.calc_dw(self.variables["w"], self.variables["u"], 
                                          self.parameters["u_c"], self.parameters["tau_w_m"], self.parameters["tau_w_p"])
        
        J_fi = ops.calc_Jfi(self.variables["u"], self.variables["v"], 
                            self.parameters["u_c"], self.parameters["tau_d"])
        J_so = ops.calc_Jso(self.variables["u"], self.parameters["u_c"],
                            self.parameters["tau_o"], self.parameters["tau_r"])
        J_si = ops.calc_Jsi(self.variables["u"], self.variables["w"],
                            self.parameters["k"], self.parameters["uc_si"], self.parameters["tau_si"])

        self.variables["u"] += self.dt * (ops.calc_rhs(J_fi, J_so, J_si) + 
                                        sum(stim.stim(i * self.dt) for stim in self.stimulations))

    def run(self, t_max: float):
        """
        Run the simulation up to time t_max.
        
        Parameters
        ----------
        t_max : float
            Maximum simulation time.
        """
        n_steps = int(round(t_max/self.dt))
        for i in range(n_steps):
            self.step(i)
            for s in self.variables:
                self.history[s].append(self.variables[s])