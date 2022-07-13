
% Example 2. Doubly Coated Disc

% Create object of type Disc, with associated folder
mydisc2 = DoublyCoatedDisc();
mydisc2 = mydisc2.set_folder('..\myfolder2');
mydisc2 = mydisc2.reinitialise();   % in case you run this script more than once

% Set properties:

% 1) Geometry properties

mydisc2.substrate.diameter = 0.0254*2;
mydisc2.substrate.thickness = .2e-3;

mydisc2.coating.thickness = 1e-6;

% 2) Material properties:

mydisc2.substrate.material = Material('silicon');
mydisc2.substrate.material.flag_anisotropic = 'single_crystal';
mydisc2.substrate.material.orientation_axis_1 = [0 0 1];	

mydisc2.coating.material = Material('silica');

% 3) Set the temperature

mydisc2.substrate.material = mydisc2.substrate.material.set_temperature(200);
mydisc2.coating.material = mydisc2.coating.material.set_temperature(200);

% 4) Mesh settings

mydisc2.mesh_settings.method = 'sweep';
mydisc2.mesh_settings.divisions = 5;  % no. of elements along the thickness                  
mydisc2.mesh_settings.fineness = mydisc2.substrate.diameter/20;  % element size   

% 5) Calculation settings

mydisc2.calculation_settings.max_number_of_modes = 10;
mydisc2.calculation_settings.mode_list = [2:5];

% 6) Calculate with Ansys

workbench_out = mydisc2.ansys_workbench_calculate();
mydisc2 = mydisc2.ansys_apdl_calculate();

% 7) Inspect Ansys results for one specific mode (Mode 4)

% Refer to the mode by the position occupied in the Ansys's list above;
% i.e. mode_number = 4 selects the 4th mode found by Ansys, 
% and stored in sub-folder "4",
% corresponding to the third mode in the vector list mode_list
% and to the third position in all vectorial results loaded into mydisc2 
% (frequency_modes; elastic_energy_modes; etc.)

mode_number = 4;   

element_volumes = mydisc2.get_element_volume(mode_number);
element_volume_changes = mydisc2.get_element_volumechange(mode_number);

[strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = ...
        mydisc2.get_element_strain(mode_number);
[stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = ...
        mydisc2.get_element_stress(mode_number);

elastic_energy = mydisc2.get_element_elastic_energy(mode_number);
[elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = ...
    mydisc2.get_element_energy_approx(mode_number);

% 8) Calculate Thermoelastic Dissipation with Vengallatore's model

mydisc2 = mydisc2.calculate_TED_vengallatore(2);
mydisc2 = mydisc2.calculate_dissipation_TED_plus_coating();

% 9) Inspect results	

disp(mydisc2.calculation_results);

% 10) Save the MatLab variable "mydisc2" (within the project folder)

mydisc2.save();




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
