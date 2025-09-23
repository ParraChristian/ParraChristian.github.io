---
title: Aeroloads
layout: default
filename: AeroLoads.md
remote_theme: pages-themes/cayman@v0.2.0
plugins:
- jekyll-remote-theme # add this line to the plugins list if you already have one
--- 
# Wing Bending, Torsion, and Shear Analysis Program

## Project Overview

This MATLAB-based analysis tool evaluates structural performance of aircraft wings under combined torsion, 
bending, and shear loads. It enables rapid assessment of wing structural integrity, providing key insights into stress
distribution, deflection, and safety margins for conceptual and preliminary design phases.

## How the Code Works

### User Inputs:

• Wing geometry (span, chord, airfoil thickness, spar/rib layout)

• Material properties (Young’s modulus, shear modulus, density)

• Applied loading conditions (lift distribution, shear forces, torsional moments, and bending moments)

### Analysis Modules:

1. Bending Analysis

   • Computes bending stress distribution using beam theory.

   • Determines maximum deflection and stress along the span.

2. Shear Analysis

   • Calculates shear flow and shear stresses in wing spars and skin panels.

   • Evaluates load paths through multi-cell wing box structures.

3. Torsional Analysis

   • Solves for twist distribution along the wing span.

   • Accounts for warping constraints and stiffness of closed vs. open sections.

### Outputs:

• Stress and deflection plots across the span.

• Shear flow distribution and torsional twist results.

• Structural safety checks based on material limits and user-defined safety factors.

### Impact

This program streamlines early-stage wing structural analysis, eliminating manual calculations and providing a 
reusable framework for rapid design iteration. It demonstrates proficiency in structural mechanics, MATLAB programming,
and aerospace system design.

## Diagram of Wing profile
<img src="/docs/AeroLoads/AeroLoads_diagram.png" width="50%"><br/>
## Sample input file
<img src="/docs/AeroLoads/AeroLoads_in.png" width="50%"><br/>
## Sample  output file
<img src="/docs/AeroLoads/AeroLoads_out1.png" width="50%"><br/>
<img src="/docs/AeroLoads/AeroLoads_out2.png" width="50%"><br/>
<img src="/docs/AeroLoads/AeroLoads_out3.png" width="50%"><br/>
<img src="/docs/AeroLoads/AeroLoads_out4.png" width="50%"><br/>
<img src="/docs/AeroLoads/AeroLoads_out5.png" width="50%"><br/>
<img src="/docs/AeroLoads/AeroLoads_out6.png" width="50%"><br/>

