## Brain Resource Localization

This project aims to address the challenge of accurately locating brain resources using a three-layer spherical head model based on Maxwell's equations.

### Description

The project involves several steps:

- **a)** Generating all possible dipoles with a resolution of one centimeter, visualizing them in 3D space, and computing the gain matrix.
- **b)** Drawing the electrode locations and labeling each electrode.
- **c)** Placing a dipole at a random location on the surface of the sphere with its direction in the radial direction.
- **d)** Simulating spiking activity (non-convulsive epilepsy) on selected bipolars, calculating the potentials in 21 electrodes, and visualizing the potentials over time for each electrode.
- **e)** Identifying the time of occurrence of positive peaks in electrode potentials, extracting windows around these peaks, and visualizing the potentials across electrodes with a color spectrum.
- **f)** Repeating part (c) using the D3_Potential_Display.m function.
- **g)** Applying MNE and WMNE algorithms to solve the inverse problem.
- **h)** Estimating dipole locations and directions for each positioning method.
- **i)** Evaluating the error in estimating dipole location and direction.
- **j)** Comparing results for surface and deep dipoles.
- **k)** Applying two selected non-parametric methods and repeating previous steps.
- **l)** Solving the inverse problem using a parametric model and a search algorithm to estimate dipole location and direction.
- **m)** Visualizing a set of adjacent dipoles and their direction vectors.
- **n)** Attributing spike activities to each bipolar and repeating relevant steps.
- **o)** Obtaining estimated moment ranges for each dipole for each positioning method.
- **p)** Drawing ROC curves for each method and comparing diagnostic accuracy.
- **q)** Repeating steps (g) to (h) with two additional non-parametric methods.

### How to Use

To utilize this project:

1. Clone the repository.
2. Run the provided scripts, ensuring necessary dependencies are installed.
3. Follow the instructions in each section of the project to perform specific tasks and analyze results.

### Files Included

- **ForwardModel_3shell.m**: Function to calculate dipole locations and the gain matrix.
- **Display_Potential_3D.m**: Function to visualize EEG potentials on the head.
- **Interictal.mat**: Matrix containing spike activities.
- **ElecPosXYZ.mat**: Matrix containing electrode positions.
- **ElecPatch.mat**: Matrix containing electrode patches.

### Dependencies

- MATLAB
- Additional MATLAB toolboxes as required.

### Contributors

- [Arya Koureshi]

### License

This project is licensed under the [License Name] License - see the [LICENSE.md](link) file for details.

### Acknowledgements

Special thanks to [Acknowledged Individuals or Organizations] for their contributions and support.

Feel free to reach out with any questions or feedback!
