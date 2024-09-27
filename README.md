# Fault Slip Potential (FSP) 2.0

### A Program for Probabilistic Estimation of Fault Slip Potential Resulting from Fluid Injection

---

### Overview

FSP 2.0 is a tool developed to screen faults near injection wells and estimate the potential for fault slip induced by fluid injection. It uses probabilistic methods to estimate the cumulative probability of a fault slipping due to increased pore pressure.

### Features

- **Probabilistic and Deterministic Models**: FSP provides both deterministic and Monte Carlo probabilistic approaches for estimating fault slip.
- **Input Parameters**: Supports input of fault strike, dip, well locations, injection rates, hydrologic parameters, and mechanical stress state parameters.
- **Hydrology Model**: Includes a simplified radial flow hydrology model, allowing users to relate injection to pore pressure changes.
- **Monte Carlo Simulations**: Uses Monte Carlo analysis to calculate probabilities of fault slip as a function of pore pressure increase.

### How It Works

1. **Mohr-Coulomb Slip Criteria**: FSP calculates pore pressure to slip on each fault using deterministic geomechanical modeling.
2. **Monte Carlo Analysis**: Probabilistic simulations are run to yield the probability of each fault slipping, considering uncertainties in input parameters.
3. **Hydrology Model**: The hydrology model assesses specific injection scenarios, providing pore pressure changes to be used in the fault slip analysis.
4. **Results**: FSP produces visual outputs including Mohr diagrams, fault maps, and cumulative distribution functions (CDF) of fault slip probability.

### Installation

FSP 2.0 requires a 64-bit Windows machine and the MATLAB Compiler Runtime (MCR) version 8.0. which is included in the [download folder](https://www.beg.utexas.edu/files/texnet-cisr/FSP%20Tutorials%20and%20Example%20Analysis.zip).

What's included in this download:

1. This README file

2. FSP_User_Guide.pdf, the code documentation

3. MatlabCompilerRuntime_windows64Bit_installer.exe, the installer for MATLAB compiler runtime: install this on your 64 bit windows machine to be able to run the FSP program. 

4. FaultSlipPotential_2.0.exe , The FSP program, double click after MCR is installed. 

5. Folder with 3 example data formatted files, for if you want to load your own faults, wells or pressure models into the FSP program, format them like these 3 files.

MCR is can also be downloaded separately from the [Mathworks website](https://www.mathworks.com/supportfiles/MCR_Runtime/R2012b/MCR_R2012b_win64_installer.exe).

1. Install the MATLAB Compiler Runtime (MCR).
2. Run the `FaultSlipPotential_2.0.exe` file to launch the application.


### Usage

A video on how to download, install, and use FSP can be found here: https://www.youtube.com/watch?v=lY0R-mt10zg

### Licensing

FSP 2.0 is licensed under the **Apache License Version 2.0**. Please refer to the `LICENSE` file for more details.

### Notes

We are currently working on a new version of FSP. Subscribe to our [newsletter](https://texnet.beg.utexas.edu/tools/fsp) to stay updated.


---



