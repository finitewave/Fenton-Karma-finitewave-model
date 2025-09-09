"""
Example script to run a 0D model simulation and plot the results.

This script sets up a simple stimulation protocol, runs the simulation,
and plots the membrane potential over time.
"""

import numpy as np
import matplotlib.pyplot as plt

from implementation.fenton_karma_0d import FentonKarma0D, Stimulation


stimulations = [Stimulation(t_start=0.1, duration=0.2, amplitude=1.0)]
t_max = 400.0

model = FentonKarma0D(dt=0.01, stimulations=stimulations)
model.run(t_max=t_max)

time = np.arange(0, t_max, model.dt)
plt.plot(time, model.history['u'])
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (u)')
plt.title('Fenton-Karma Model Simulation')
plt.grid()
plt.show()

