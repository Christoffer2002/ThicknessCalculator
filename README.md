#Swanepoel optical thickness and constants calculator
This project is an implementation of the Swanepoel method utilizing the transmission spectrum through a thin film to determine the thickness and refractive index of said film. This script follows the R. Swanepoel approach descriped in ** R Swanepoel 1983 J. Phys. E: Sci. Instrum. 16 1214**.

#Method overview
The Swanepoel method utilizes the positive and negative intereference fringes within the transmission spectrum to calculate film thickness.
The key idea is to extract two smooth envelopes, one that follows the maxima T_m(λ) and one following the minima T_m(\lambda).
Several different strategies exist for constructing these envelopes. In this pipeline, maxima and minima are first detected using a peak-finding algorithm, after which the extracted peak sequences are smoothed using a low-order polynomial fit.
This produces a pair of continuous envelope curves that trace the overall upper and lower boundaries of the fringe pattern while suppressing noise and local fluctuations.
Once 


#Required Input
The pipeline expects a transmission spectrum with:
* Wavelength (in nm)
* Transmission (0-100)

#Physical setup


#Installation
```
git clone https://github.com/Christoffer2002/ThicknessCalculator.git
cd swanepoel-calculator
pip install -r requirements.txt
```
