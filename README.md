Binary Black Hole Inspiral Simulation

Overview

This MATLAB script simulates the inspiral of a binary black hole system, incorporating:

Post-Newtonian corrections up to the 2.5PN order

Gravitational wave energy loss calculations

Numerical integration of orbital dynamics

Generation of a gravitational waveform

Features

Accurate Orbit Modeling: Uses Newtonian and post-Newtonian corrections to model black hole inspirals.

Gravitational Wave Generation: Computes the emitted gravitational wave strain using the quadrupole formula.

Visualization: Produces plots of the binary orbit and corresponding gravitational waveforms.

Requirements

MATLAB

Basic knowledge of orbital dynamics and relativity (optional)

How to Run

Open MATLAB and navigate to the directory containing binary_black_hole.m.

Run the script using:

binary_black_hole

The simulation will generate:

A plot of the inspiraling orbit

A plot of the gravitational waveform

Functions

equations: Computes the equations of motion with post-Newtonian corrections.

compute_waveform: Generates the gravitational wave strain from the binary system.

Notes

The simulation assumes a simplified 2-body problem with no external perturbations.

The merger phase is approximated and does not include full numerical relativity corrections.

Future Improvements

Implement higher-order PN terms for better accuracy.

Include a waveform model that captures the ringdown phase.

Introduce real-time animation of the inspiral process.

References

Blanchet, L. (2014). Gravitational Radiation from Post-Newtonian Sources.

Misner, C., Thorne, K., & Wheeler, J. (1973). Gravitation.
