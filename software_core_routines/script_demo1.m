
% Example 1. Uncoated Disc

% This demo uses options 2-b, 3-a, 4-a from the documentation.

% Create object of type Disc, with associated folder
mydisc1 = Disc();
mydisc1 = mydisc1.set_folder('..\myfolder1');
mydisc1 = mydisc1.reinitialise();   % in case you run this script more than once

% Set properties:

% 1) Geometry properties (thickness = 300 micron and diameter = 1 inch)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
mydisc1.substrate.diameter = 0.0254;
mydisc1.substrate.thickness = .3e-3;

% 2) Material properties

mydisc1.substrate.material = Material('silicon');

mydisc1.substrate.material.flag_anisotropic = 'single_crystal';
mydisc1.substrate.material.orientation_axis_1 = [1 1 1];	

% 3) Set the temperature

mydisc1.substrate.material = mydisc1.substrate.material.set_temperature(300);

% 4) Mesh settings

mydisc1.mesh_settings.method = 'sweep';
mydisc1.mesh_settings.divisions = 5;  % no. of elements along the thickness     

mydisc1.mesh_settings.fineness = mydisc1.substrate.diameter/20;  % element size   

% 5) Calculation settings

% Let Ansys create a list of modes, according to the following constraints
mydisc1.calculation_settings.max_number_of_modes = 20;
mydisc1.calculation_settings.max_frequency = 100000;
mydisc1.calculation_settings.min_frequency = 100;  % already the default	

% Create a vector-list of the modes of interest, whose results will be saved.
% (The numbers in this vector-list are referred to the position occupied in the 
% Ansys's list above; i.e. mode_list = [2,4] selects the 2nd and the 4th mode 
% found by Ansys, within the constraints specified)

mydisc1.calculation_settings.mode_list = [2,4];  % modes of interest

% The results for each mode will be saved in a sub-folder labeled by 
% the mode number, e.g. "2" and "4" in this case

% 6) Calculate with Ansys

workbench_out = mydisc1.ansys_workbench_calculate();  % prepare input file
mydisc1 = mydisc1.ansys_apdl_calculate();  % edit input file and execute

% 7) Inspect Ansys results for one specific mode (Mode 4)

% Refer to the mode by the position occupied in the Ansys's list above;
% i.e. mode_number = 4 selects the 4th mode found by Ansys, 
% and stored in sub-folder "4",
% corresponding to the second mode in the vector list mode_list
% and to the second position in all vectorial results loaded into mydisc1 
% (frequency_modes; elastic_energy_modes; etc.)

mode_number = 4;   

element_volumes = mydisc1.get_element_volume(mode_number);
element_volume_changes = mydisc1.get_element_volumechange(mode_number);

[strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = ...
        mydisc1.get_element_strain(mode_number);
[stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = ...
        mydisc1.get_element_stress(mode_number);

elastic_energy = mydisc1.get_element_elastic_energy(mode_number);
[elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = ...
    mydisc1.get_element_energy_approx(mode_number);

% 8) Calculate Thermoelastic Dissipation according to Vengallatore's model

mydisc1 = mydisc1.calculate_TED_vengallatore();
mydisc1 = mydisc1.calculate_dissipation_TED();

% 9) Inspect results	

disp(mydisc1.calculation_results);

% 10) Save the MatLab variable "mydisc1" (within the project folder)

mydisc1.save();




% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

% F. Fabrizi (2022); federica.fabrizi@uniurb.it





