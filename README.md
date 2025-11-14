# Swanepoel optical thickness and constants calculator
This project is an implementation of the Swanepoel method utilizing the transmission spectrum through a thin film to determine the thickness and refractive index of said film. This script follows the R. Swanepoel approach descriped in ** R Swanepoel 1983 J. Phys. E: Sci. Instrum. 16 1214**.

# How to use
The ```examples/``` folder contains simple test cases demonstrating how to use the pipeline with real data. ```examples/run_single.py``` is an example of a script running the full pipeline on one single optical measurement and ```examples/run_multiple.py``` is an example of multiple measurements.

# Method overview
The Swanepoel method utilizes the positive and negative intereference fringes within the transmission spectrum to calculate film thickness.
The key idea is to extract two smooth envelopes, one that follows the maxima $T_M(\lambda)$ and one following the minima $T_m(\lambda)$.
Several different strategies exist for constructing these envelopes. In this pipeline, maxima and minima are first detected using a peak-finding algorithm, after which the extracted peak sequences are smoothed using a low-order polynomial fit.
This produces a pair of continuous envelope curves that trace the overall upper and lower boundaries of the fringe pattern while suppressing noise and local fluctuations.
Once the $T_m(\lambda)$ and $T_M(\lambda)$ envelopes are built a range of equations can be applied allowing computation of the wavelength-dependent refractive index $n(\lambda)$ and the film thickness.


# Required Input
The pipeline is built to accept csv files from the F20 Thin film analyzer and therefor expects a transmission spectrum with:
* Wavelength (in nm)
* Transmission (in % -> 0-100)

# Physical setup
The Swanepoel method depends critically on the optical coherence conditions of the measurement setup. While the mathematics assumes perfectly coherent reflection paths, real light sources have finite coherence lengths, which impose practical constraints. The next subchapters introduce possible constraints in the physical setup that can introduce high or low frequency noise into the transmission spectrum and result in wrong thickness predictions. 

## Coherence Length and Fringe Visibility
The coherence length describes the distance over which the light maintains a stable phase relationship. The coherence length must be longer than the optical path difference created by the film.
Otherwise, reflections from deeper layers lose phase correlation, and the measured transmission becomes smooth.
For interference fringes to form:
## Film Thickness Constraint
If the coherence length is shorter than the film thickness, interference patterns disappear entirely, making the Swanepoel method unusable.
A fringe-free spectrum does not carry enough information to extract film thickness.
## Substrate Constraint
The coherence length must also be shorter than the substrate thickness.
If it exceeds the substrate thickness, the substrate itself forms a Fabry–Perot cavity and produces its own fringes. These substrate fringes overlap the film fringes, making it impossible to isolate the film’s contribution reliably.
## Wavelength Dependence
Coherence length increases with wavelength. At longer wavelengths it may eventually exceed substrate thickness or create high-frequency oscillations. These effects can distort envelope extraction and introduce incorrect fringe order assignment.
In summary, the coherence length must lie between the film thickness and the substrate thickness:
$d<L_{coherence}<d_{substrate}$
## Tip 👍
Use ```frequency.py``` after an initial thickness estimate to verify that your spectral sampling meets Nyquist for the predicted fringe rate; if not, increase spectral resolution or restrict the analysis range
# Installation
```
git clone https://github.com/Christoffer2002/ThicknessCalculator.git
cd ThicknessCalculator
pip install .
```
