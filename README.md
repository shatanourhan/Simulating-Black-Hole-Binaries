# Binary Black Hole Inspiral Simulation

## Overview

This MATLAB script simulates the inspiral of a binary black hole system, incorporating:

- Post-Newtonian corrections up to the 2.5PN order
- Gravitational wave energy loss calculations
- Numerical integration of orbital dynamics
- Generation of a gravitational waveform

## Features

- **Accurate Orbit Modeling**: Uses Newtonian and post-Newtonian corrections to model black hole inspirals.
- **Gravitational Wave Generation**: Computes the emitted gravitational wave strain using the quadrupole formula.
- **Visualization**: Produces plots of the binary orbit and corresponding gravitational waveforms.

## Requirements

- MATLAB
- Basic knowledge of orbital dynamics and relativity (optional)

3. The simulation will generate:
   - A plot of the inspiraling orbit
   - A plot of the gravitational waveform

## Functions

- `equations`: Computes the equations of motion with post-Newtonian corrections.
- `compute_waveform`: Generates the gravitational wave strain from the binary system.

## Equations Used

### 1. Newtonian Equations of Motion:
$$
\frac{d\mathbf{r}}{dt} = \mathbf{v}
$$
$$
\frac{d\mathbf{v}}{dt} = - \frac{G(m_1 + m_2)}{r^3} \mathbf{r}
$$

### 2. Post-Newtonian Corrections (Up to 2.5PN Order):
$$
\frac{d\mathbf{v}}{dt} = -\frac{G(m_1 + m_2)}{r^3} \mathbf{r} + \text{PN corrections}
$$

### 3. Gravitational Wave Energy Loss (Peters & Mathews Formula):
$$
\frac{dE}{dt} = -\frac{32}{5} \frac{G^4}{c^5} \frac{m_1^2 m_2^2 (m_1 + m_2)}{r^5}
$$

### 4. Gravitational Wave Strain (Quadrupole Approximation):
$$
 h_{+}(t) \propto \frac{2G}{c^4} \frac{I_{xx} - I_{yy}}{D}
$$

where \( I_{xx} \) and \( I_{yy} \) are components of the mass quadrupole tensor, and \( D \) is the distance to the observer.

## Notes

- The simulation assumes a simplified 2-body problem with no external perturbations.
- The merger phase is approximated and does not include full numerical relativity corrections.



## References

- Blanchet, L. (2014). Gravitational Radiation from Post-Newtonian Sources.
- Misner, C., Thorne, K., & Wheeler, J. (1973). Gravitation.
