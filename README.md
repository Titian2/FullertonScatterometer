# Fullerton Scatterometer Analysis Framework

This directory contains MATLAB scripts and functions that build on the framework of the Fullerton Scatterometer to analyze scattering data and estimate the Bidirectional Reflectance Distribution Function (BRDF) from an in-situ scattering system. These codes are adapted from earlier work by Joshua Smith at UC Fullerton, specifically the "Csuflog20140820 toy model."

Link to Fullerton Log book : https://wiki.ligo.org/OPT/Csuflog20140820


## Overview of Functions

### `Model_Play.m`
This script implements several toy models to demonstrate why the Total Integrated Scatter (TIS) estimation works. It includes:
- Circular ROI-based models for analyzing scattering data.
- Dynamic simulations of scattering systems with evolving background brightness and growing circles.
- Visualization of ROI placement and live updates of fit parameters.

Key features:
- Linear fitting of sum vs. area for circular ROIs.
- Dynamic updates of scattering data with live plots of y-intercept evolution.
- Nested circular ROI analysis for complex scattering systems.

### `interactive_circle.m`
This function provides an interactive tool for visualizing and adjusting circular ROIs on an image. Key features include:
- Loading and displaying an image.
- Interactive resizing of circular ROIs using mouse and keyboard inputs.
- Visualization of the ROI placement in real time.

### `runDynamicTIS.m`
This function simulates a dynamic ROI-based scattering system with arbitrary parameters. Key features include:
- Randomized placement of circles with adjustable parameters such as size, amplitude distribution, and placement bias.
- Generation of output TSV files summarizing the scattering data.
- Reproducible results through seeding.

### `StudyPlotting_and_Clustering.m`
This script provides tools for analyzing and visualizing the results of multiple scattering experiments. Key features include:
- Plotting y-intercept curves grouped by experiment class, with color-coded outcomes (e.g., negative or positive intercepts).
- Scatter plots of final y-intercepts by experiment class, with jittered points for better visualization.
- Summary statistics for each experiment class, such as the proportion of negative outcomes.

### `TIS_Controller.m`
This script serves as the main controller for running and managing TIS-related simulations and analyses. Key features include:
- Configuring and executing multiple scattering simulations with varying parameters.
- Managing input and output files for experiments.
- Automating the generation of summary data for further analysis.

## Usage

1. **Run `Model_Play.m`** to explore toy models and visualize the relationship between ROI area and scattering data.
2. **Use `interactive_circle.m`** to interactively define and adjust circular ROIs on your scattering images.
3. **Call `runDynamicTIS.m`** to simulate scattering systems with custom parameters and generate summary data for further analysis.
4. **Use `StudyPlotting_and_Clustering.m`** to analyze and visualize the results of multiple experiments.
5. **Run `TIS_Controller.m`** to automate and manage TIS simulations and analyses.

## Dependencies

- MATLAB R2020b or later.
- Image Processing Toolbox (for `interactive_circle.m`).

## Acknowledgments

These codes are based on the foundational work by Joshua Smith at UC Fullerton and have been adapted to extend the analysis methods for the Fullerton Scatterometer.

## Future Work

This framework will continue to evolve to improve BRDF estimation methods and integrate additional scattering analysis techniques.