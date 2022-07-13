
% Example 4. Cantilever Fiber

% Create object of type CantileverFiber, with associated folder

myfiber = CantileverFiber();
myfiber = myfiber.set_folder('..\myfolder4');
myfiber = myfiber.reinitialise();   % in case you run this script more than once

% Set properties:

% 1) Geometry properties:

myfiber.substrate.diameter = 2e-3;     % [m] 2 mm
myfiber.substrate.thickness = 335e-3;  % [m] 335 mm

% 2) Material properties:

myfiber.substrate.material = Material('silica');

% 3) Set the temperature:

myfiber.substrate.material = myfiber.substrate.material.set_temperature(300);

% 4) Mesh settings

myfiber.mesh_settings.method = 'non-sweep';
myfiber.mesh_settings.fineness = 0;  % elready the default   

% 5) Calculation settings

% Let Ansys create a list of modes, according to the following constraints
myfiber.calculation_settings.max_number_of_modes = 20;
myfiber.calculation_settings.max_frequency = 100000;
myfiber.calculation_settings.min_frequency = 100;  % already the default

% Create a vector-list of the modes of interest, whose results will be saved.
% (The numbers in this vector-list are referred to the position occupied in the 
% Ansys's list above; i.e. mode_list = [2,4] selects the 2nd and the 4th mode 
% found by Ansys, within the constraints specified)

myfiber.calculation_settings.mode_list = [2,4];  % modes of interest

% The results for each mode will be saved in a sub-folder labeled by 
% the mode number, e.g. "2" and "4" in this case

% 6) Calculate with Ansys

workbench_out = myfiber.ansys_workbench_calculate();
myfiber = myfiber.ansys_apdl_calculate();

% 8) Calculate Thermoelastic Dissipation according to Zener's model

myfiber = myfiber.calculate_TED_zener();
myfiber = myfiber.calculate_dissipation_TED();

% 9) Inspect results	

disp(myfiber.calculation_results);

% 10) Save the MatLab variable "myfiber" (within the project folder)

myfiber.save();




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

