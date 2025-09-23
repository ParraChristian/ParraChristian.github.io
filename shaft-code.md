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

<summary> Project MATLAB Code <summary>

### You can add a header

```matlab
   % MAE 190 
% Steel Shaft Design Optimization Tool
% Christian Parra
% A16402854
% 3/17/2024
% all calculations are made in imperial units[lb,in,ksi/psi]
%-------------------------------------------------------------------------

% user interface for shaft criteria
function mae190test()

 comp=uifigure('Name',"Steel Shaft Design Optimization Tool");
 comp.Position = [1 1 905 756];

 % Create GridLayout
 g = uigridlayout(comp);
 g.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
 g.RowHeight = {25, 40, '1x', '1x', 40, '1x', '1x', '1x'};

 % Create LabeloftextareaTextAreaLabel
 ga = uilabel(g);
 ga.HorizontalAlignment = 'center';
 ga.Layout.Row = 1;
 ga.Layout.Column = 1;
 ga.Text = 'Christian Parra';

 % Create MAE190APPTextArea
 gb = uitextarea(g);
 gb.HorizontalAlignment = 'center';
 gb.Placeholder = 'text area place holder';
 gb.Layout.Row = 1;
 gb.Layout.Column = [2 6];
 gb.Value = {'Program to calculate the minimum diamiter for a steel shaft with a shoulder fillet'};

 % Create MaterialSelectionLabel
 gc = uilabel(g);
 gc.HorizontalAlignment = 'center';
 gc.Layout.Row = 2;
 gc.Layout.Column = [1 3];
 gc.Text = 'Material Selection';

 % Create TabGroup2
 gd = uitabgroup(g);
 gd.Layout.Row = [3 4];
 gd.Layout.Column = [1 3];

 % Create CustomTab
 ge = uitab(gd);
 ge.Title = 'Custom';

 % Create TensileStrengthKsiEditFieldLabel
 gf = uilabel(ge);
 gf.HorizontalAlignment = 'right';
 gf.Position = [24 134 91 30];
 gf.Text = {'Tensile Strength'; '[Ksi]'};

 % Create TensileStrengthKsiEditField
 gg = uieditfield(ge, 'numeric');
 gg.Limits = [50 250];
 gg.Position = [130 142 100 22];
 gg.Value = 50;

 % Create YieldStrengthKsiEditFieldLabel
 gh = uilabel(ge);
 gh.HorizontalAlignment = 'right';
 gh.Position = [35 76 80 30];
 gh.Text = {'Yield Strength'; '[Ksi]'};

 % Create YieldStrengthKsiEditField
 gi = uieditfield(ge, 'numeric');
 gi.Limits = [50 250];
 gi.Position = [130 84 100 22];
 gi.Value = 50;

 % Create SurfaceFinishListBox_2Label
 gj = uilabel(ge);
 gj.HorizontalAlignment = 'right';
 gj.Position = [288 166 81 22];
 gj.Text = 'Surface Finish';

 % Create SurfaceFinishListBox_2
 gk = uilistbox(ge);
 gk.Position = [282 76 150 90];
 gk.Items={'Ground','Machined or Cold drawn','Hot-rolled','As-Forged'};

 % Create DatabaseTab
 gl = uitab(gd);
 gl.Title = 'Database';

 % Create MaterialSelectionListBoxLabel
 gm = uilabel(gl);
 gm.HorizontalAlignment = 'right';
 gm.Position = [53 169 101 22];
 gm.Text = 'Material Selection';

 % Create MaterialSelectionListBox
 gn = uilistbox(gl);
 gn.Items = {'1006HR', '1006CD', '1010HR', '1010CD', '1015HR', '1015CD', '1018HR', '1018CD', '1020HR', '1020CD', '1030HR', '1030CD', '1035HR', '1035CD', '1040HR', '1040CD', '1045HR', '1045CD', '1050HR', '1050CD1', '060HR', '1080HR', '1095HR'};
 gn.Position = [169 54 140 139];
 gn.Value = '1006HR';

 % Create DropDown
 go = uidropdown(g);
 go.Items = {'Custom', 'Database'};
 go.Layout.Row = 2;
 go.Layout.Column = 1;
 go.Value = 'Custom';

 % Create LoadslbinPanel
 gp = uipanel(g);
 gp.TitlePosition = 'centertop';
 gp.Title = 'Loads [lb-in]';
 gp.Layout.Row = [5 6];
 gp.Layout.Column = [1 3];

 % Create MaxMomentEditFieldLabel
 gq = uilabel(gp);
 gq.HorizontalAlignment = 'right';
 gq.Position = [13 86 76 22];
 gq.Text = 'Max Moment';

 % Create MaxMomentEditField
 gr = uieditfield(gp, 'numeric');
 gr.Position = [104 86 100 22];
 gr.Limits = [-60000 60000];

 % Create MinMomentEditFieldLabel
 gs = uilabel(gp);
 gs.HorizontalAlignment = 'right';
 gs.Position = [228 86 73 22];
 gs.Text = 'Min Moment';

 % Create MinMomentEditField
 gt = uieditfield(gp, 'numeric');
 gt.Position = [316 86 100 22];
 gt.Limits = [-60000 60000];

 % Create MaxTorqueEditFieldLabel
 gu = uilabel(gp);
 gu.HorizontalAlignment = 'right';
 gu.Position = [21 38 68 22];
 gu.Text = 'Max Torque';

 % Create MaxTorqueEditField
 gv = uieditfield(gp, 'numeric');
 gv.Position = [104 38 100 22];
 gv.Limits = [-60000 60000];

 % Create MinTorqueEditFieldLabel
 gw = uilabel(gp);
 gw.HorizontalAlignment = 'right';
 gw.Position = [236 38 65 22];
 gw.Text = 'Min Torque';

 % Create MinTorqueEditField
 gx = uieditfield(gp, 'numeric');
 gx.Position = [316 38 100 22];
 gx.Limits = [-60000 60000];

 % Create OptionsPanel
 gy = uipanel(g);
 gy.TitlePosition = 'centertop';
 gy.Title = 'Options';
 gy.Layout.Row = [7 8];
 gy.Layout.Column = [1 3];

 % Create ReliabilityListBoxLabel
 gz = uilabel(gy);
 gz.HorizontalAlignment = 'right';
 gz.Position = [39 170 56 22];
 gz.Text = 'Reliability';

 % Create ReliabilityListBox
 gaa = uilistbox(gy);
 gaa.Items = {'50%', '90%', '95%', '99%', '99.9%', '99.99%'};
 gaa.Position = [30 20 74 138];
 gaa.Value = '50%';

 % Create SafetyFactorEditFieldLabel
 gbb = uilabel(gy);
 gbb.HorizontalAlignment = 'right';
 gbb.Position = [283 59 77 22];
 gbb.Text = 'Safety Factor';

 % Create SafetyFactorEditField
 gcc = uieditfield(gy, 'numeric');
 gcc.Limits = [1 Inf];
 gcc.HorizontalAlignment = 'center';
 gcc.Position = [272 28 100 22];
 gcc.Value = 1;

 % Create rdListBoxLabel
 gdd = uilabel(gy);
 gdd.HorizontalAlignment = 'center';
 gdd.Position = [313 170 25 22];
 gdd.Text = 'r/d';

 % Create rdListBox
 gee = uilistbox(gy);
 gee.Items = {'.02(Sharp)', '.1(Well Rounded)'};
 gee.Position = [245 103 161 55];
 gee.Value = '.02(Sharp)';

 % Create DdListBoxLabel
 gff = uilabel(gy);
 gff.HorizontalAlignment = 'right';
 gff.Position = [163 170 25 22];
 gff.Text = 'D/d';

 % Create DdListBox
 ggg = uilistbox(gy);
 ggg.Items = {'1.10', '1.15', '1.20', '1.25', '1.30'};
 ggg.Position = [137 20 77 138];
 ggg.Value = '1.10';

 % Create ResultsPanel
 ghh = uipanel(g);
 ghh.TitlePosition = 'centertop';
 ghh.Title = 'Results';
 ghh.Layout.Row = [3 6];
 ghh.Layout.Column = [4 6];

 % Create RunCalculationsButton
 gii = uibutton(g, 'push');
 gii.ButtonPushedFcn = matlab.apps.createCallbackFcn(comp, ...
     @RunCalculationsButtonPushed, true);
 gii.Layout.Row = [7 8];
 gii.Layout.Column = [4 6];
 gii.Text = 'Run Calculations';
 RESULT={''};
 gjj = uilabel(ghh, 'Position', [100 100 500 200],'Text',RESULT);



    function RunCalculationsButtonPushed(~,~)

        %custom vs dictionary input
        k=go.ValueIndex;

        %select sy and su based on custom vs database
        if k==1
            Su = gg.Value;
            Sy = gi.Value;

            %Material Finish Table
            Finish_a=[1.34 2.7 14.4 39.9];
            Finish_b=[-.085 -.265 -.718 -.955];

            j=gk.ValueIndex;
            ka=Finish_a(j)*Su^Finish_b(j);

        else
            %material dictionary
            Material_names=['1006HR';'1006CD';'1010HR';'1010CD';'1015HR';...
                '1015CD';'1018HR';'1018CD';'1020HR';'1020CD';'1030HR';...
                '1030CD';'1035HR';'1035CD';'1040HR';'1040CD';'1045HR';...
                '1045CD';'1050HR';'1050CD';'1060HR';'1080HR';'1095HR'];
            
            Material_Su=[43; 48; 47; 53; 50; 56; 58; 64; 55; 68; 68; 76;...
                72; 80; 76; 85; 82; 91; 90; 100; 98; 112; 120];
            
            Material_Sy=[24; 41; 26; 44; 27.5; 47; 32; 54; 30; 57; 37.5;...
                64; 39.5; 67; 42; 71; 45; 77; 49.5; 84; 54; 61.5; 66];
            
            Material_table=table(Material_names,Material_Su,Material_Sy);


            indx1=gn.ValueIndex;
            Su=Material_Su(indx1);
            Sy=Material_Sy(indx1);
            ka=1;
        end

        %Size factor Kb, initial=.9
        kb=1;
        %Kc=1 Kd=1

        %reliability factor Ke
        Realiability_choice=gaa.ValueIndex;

        if Realiability_choice==1%reliability == 50
            ke=1;
        elseif Realiability_choice==2%reliability == 90
            ke=.897;
        elseif Realiability_choice==3%reliability == 95
            ke=.868;
        elseif Realiability_choice==4%reliability == 99
            ke=.814;
        elseif Realiability_choice==5%reliability == 99.9
            ke=.753;
        elseif Realiability_choice==6%reliability == 99.99
            ke=.702;
        end

        %overall endurance limit
        if Su<=200
            Se_prime=.5*Su*1000;%psi
        else
            Se_prime=100*1000;%psi
        end

        %endurance limit Se at filet initial
        Se = ka*kb*ke*Se_prime;

        %first guess Kf and Kfs
        shoulder_type=gee.ValueIndex;
        if shoulder_type == 1
            Kf=2.7; Kfs=2.2; r_d=.02;
        elseif shoulder_type == 2
            Kf=1.7; Kfs=1.5; r_d=.1;
        end

        %factor of safety
        n = gcc.Value;
        
        %impory T and M
        T_max = gv.Value;
        T_min = gx.Value;
        M_max = gr.Value;
        M_min = gt.Value;
        
        %check T and M delta is positive

         if M_min > M_max
             
             uialert(comp,"Ensure that Minimum moment is less than or equal to Maximum moment","Moment Error");
         end

         if T_min > T_max
             
             uialert(comp,"Ensure that Minimum torque is less than or equal to Maximum torque","Torque Error");
         end

        %mean and alternating T and M
        Ta=(T_max-T_min)/2;%lb-in
        Tm=(T_max+T_min)/2;%lb-in

        Ma=(M_max-M_min)/2;%lb-in
        Mm=(M_max+M_min)/2;%lb-in

        %getting first estimation of diamiter using goodman
        di(1)=(16*n/pi* (sqrt(4*(Kf*Ma)^2+3*(Kfs*Ta)^2)/Se +...
            sqrt(4*(Kf*Mm)^2+3*(Kfs*Tm)^2)/Su/1000))^(1/3);


        %now can get to iteration stage
        D_d=str2double(ggg.Value);
        %find Kt, using interpolation from fesired D/d
        if  D_d>1.2
            Kt_b=.93836 * r_d ^ -.25759;%line at 1.5
            Kt_a=.97098 * r_d ^ -.21796;%line at 1.2
            Kt=Kt_a + (Kt_b-Kt_a)/(1.5-1.2)*(D_d-1.2);
        elseif D_d<=1.2
            Kt_b=.97098 * r_d ^ -.21796;%line at 1.2
            Kt_a=.95120 * r_d ^ -.23757;%line at 1.1
            Kt=Kt_a + (Kt_b-Kt_a)/(1.2-1.1)*(D_d-1.1);
        end

        %find Kts, using interpolation from desired D/d
        if D_d>1.2
            Kts_b=.84897 * r_d ^ -.23161;%line at 1.33
            Kts_a=.83425 * r_d ^ -.21649;%line at 1.2
            Kts=Kts_a + (Kts_b-Kts_a)/(1.33-1.2)*(D_d-1.2);
        elseif D_d<=1.2
            Kts_b=.83425 * r_d ^ -.21649;%line at 1.2
            Kts_a=.90337 * r_d ^ -.12692;%line at 1.09
            Kts=Kts_a + (Kts_b-Kts_a)/(1.2-1.09)*(D_d-1.09);
        end
        %loop begins only recalculating factor that depend on the diamiter of the
        %shaft
        k=2;
        dd=1;
        %only stops after the % error between the last 2 diamiter iterations are
        %less than 0.1%

        while dd>=.001
            %find q
            a_sqrt=.246-3.08*(10^-3)*Su+1.51*(10^-5)*Su^2-2.67*(10^-8)*Su^3;
            r=r_d*di(k-1);
            q=1/(1+a_sqrt/sqrt(r));

            %find qs
            as_sqrt=.190-2.51*(10^-3)*Su+1.35*(10^-5)*Su^2-2.67*(10^-8)*Su^3;
            qs=1/(1+as_sqrt/sqrt(r));

            %finding Kf and Kfs
            Kf=1+q*(Kt-1);
            Kfs=1+qs*(Kts-1);

            %finding size factor kb accurately
            if di(k-1)<.3
                kb=(.3/.3)^-.107;
            elseif di(k-1)>=.3 && di(k-1)<=2 %in
                kb=(di(k-1)/.3)^-.107;
            elseif di(k-1)>2 && di(k-1)<=10
                kb=.92*di(k-1)^-.157;
            elseif di(k-1)>10
                kb=.92*di(k-1)^-.157;
            end

            %new endurance limit at critical point
            Se = ka*kb*ke*Se_prime;
            di(k)=(16*n/pi* (sqrt(4*(Kf*Ma)^2+3*(Kfs*Ta)^2)/Se +...
                sqrt(4*(Kf*Mm)^2+3*(Kfs*Tm)^2)/Su/1000))^(1/3);

            dd=abs((di(k)-di(k-1))/di(k-1));

            k=k+1;

            %check von mises for yielding
            Sigma_von_a=sqrt((32*Kf*Ma/pi/di(k-1)^3)^2+3*(16*Kfs*Ta/pi/di(k-1)^3));
            Sigma_von_m=sqrt((32*Kf*Mm/pi/di(k-1)^3)^2+3*(16*Kfs*Tm/pi/di(k-1)^3));
            Sigma_von_max=Sigma_von_a+Sigma_von_m;
            ny=Sy/Sigma_von_max*1000;
        end
        if ny<1
            %if yields give recomendations and ask for new inputs
            uialert(comp,"input conditions will lead to a shaft that yields. Please choose do one of the following" + ...
                " choose a different material with a higher yielding stress" + ...
                "decrease the loads that the shaft experience","Yield Error");
            gjj.Text ={''};


        else
            %if doues not yield, output athe converged diamiters
            d=di(k-1);
            D=D_d*d;

            gjj.Text = {'minimum diameter for smaller side of rod [d]';'';
                ['                            d = ',num2str(d),' in']; '';'';'';...
                'minimum diameter for larger side of rod [D]';''; ...
                ['                            D = ',num2str(D),' in']};

        
        end

    end
end
```

</details>

