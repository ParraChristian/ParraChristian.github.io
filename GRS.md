---
title: Grid Refinement Study - Vortex Cooled Rocket Engine 
layout: default
filename: GRS.md
remote_theme: pages-themes/cayman@v0.2.0
plugins:
- jekyll-remote-theme # add this line to the plugins list if you already have one
--- 
# Grid Refinement Study - Vortex Cooled Rocket Engine 

## Project Overview

This project applied a systematic grid refinement study in ANSYS Fluent to evaluate solution accuracy and establish grid independence.
By refining the mesh and analyzing changes in simulation results, I quantified numerical error, verified solver reliability, 
and ensured confidence in the CFD predictions.

## Problem Setup

Defined geometry and boundary conditions for the flow case.

Established baseline mesh, then created medium and fine meshes with consistent refinement ratios.

## Simulation in ANSYS Fluent

Ran CFD simulations across all mesh levels under the same flow conditions.

Extracted engineering quantities of interest such as pressure distribution, velocity field, and wall shear stress.

## Grid Convergence Analysis

Compared results between coarse, medium, and fine meshes.

Applied Richardson Extrapolation to estimate the asymptotic (true) solution.

Calculated Grid Convergence Index (GCI) to quantify numerical uncertainty.

## Key Findings

Verified that results converged consistently with mesh refinement.

Demonstrated grid independence for the chosen output variables.

Highlighted trade-offs between computational cost and accuracy.

## Grid Refinement Study
<img src="/docs/GRS/GRS_Study1.png" width="50%"><br/>
<img src="/docs/GRS/GRS_Study2.png" width="50%">

## Simulation Results
<img src="/docs/GRS/GRS_Results.png" width="50%">


## [Grid Refinement Study Report](docs/GRS/GRS Final Report.pdf)
