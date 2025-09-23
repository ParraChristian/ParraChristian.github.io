---
title: Shaft Code
layout: default
filename: shaft-code.md
remote_theme: pages-themes/cayman@v0.2.0
plugins:
- jekyll-remote-theme # add this line to the plugins list if you already have one
--- 
# Steel Shaft Design Optimization Tool

This MATLAB-based design tool optimizes shaft dimensions according to ASME fatigue criteria and user-defined safety factors, streamlining a traditionally time-consuming design process.<br/>

- Automated Sizing: Calculates the minimum required shaft diameter from user inputs such as torque, bending moments, alternating loads, material properties, and safety factor.<br/>

- Built-in Reliability: Incorporates automated stress concentration factor selection and ASME fatigue failure analysis, eliminating manual chart lookups while ensuring accuracy.<br/>

- Impact: Demonstrates practical application of mechanical design principles, combining coding, fatigue analysis, and optimization to deliver efficient and reliable shaft designs.<br/>


## Inputs

Material Properties<br/>
- Choose from a built-in list of ASME steel grades, or<br/>
- Manually input key properties: yield strength, tensile strength, and surface finish.<br/>

Loading Conditions<br/>
- Define the minimum and maximum moment and torque values the shaft will experience.<br/>

Additional Options<br/>
- Customize design parameters such as:<br/>
 -  Shaft reliability<br/>
 -  Factor of safety<br/>
 -  Diameter ratio (D/d)<br/>
 -  Shoulder fillet (r/d)<br/>

## Outputs
Minimum rod diamiter required to sustain applied loads for both the small and large portion of the shaft


## Default Configuration<br/>
<img width="400" height="400" alt="Image" src="https://github.com/user-attachments/assets/89e5cd04-3d8e-4645-b3be-4d729a0837ab" /><br/>
## Test Case <br/>
<img width="400" height="400" alt="Image" src="https://github.com/user-attachments/assets/4eb69f5d-1900-4947-9373-18b05c8008f0" /><br/>

<details>

<summary>Tips for collapsed sections</summary>

### You can add a header

You can add text within a collapsed section.

You can add an image or a code block, too.

```ruby
   puts "Hello World"
```

</details>

