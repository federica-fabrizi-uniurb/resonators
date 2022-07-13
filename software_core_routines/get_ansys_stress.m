
function [stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = get_ansys_stress(obj,mode_number)

    % Loads from the .txt files stored in output_files, for a given mode, the
    % element-by-element results for the stresses, as calculated by Ansys.
    % The Cartesian reference frame in which they are expressed is dependent on the element type, 
    % and may vary from one element to another (e.g. from one element belonging to the substrate
    % to one belonging to the coating), 
    % however it is guaranteed to coincide
    % with the frame for the strains on the SAME element,
    % so that the scalar energy quantities can be readily evaluated.
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    % mode_number - integer. 
    % stressXX, stressYY, etc. - arrays Nx1 double (N = number of elements).
    
    OutputFolder = obj.calculation_settings.apdl_output_folder;
    ModeFolder = fullfile(OutputFolder,num2str(mode_number));
    
    DataFile = fullfile(ModeFolder,'data_stressXX.txt');   
    stressXX = dlmread(DataFile);
    
    DataFile = fullfile(ModeFolder,'data_stressYY.txt');   
    stressYY = dlmread(DataFile);  
    
    DataFile = fullfile(ModeFolder,'data_stressZZ.txt');   
    stressZZ = dlmread(DataFile);
    
    DataFile = fullfile(ModeFolder,'data_stressXY.txt');   
    stressXY = dlmread(DataFile);    
    
    DataFile = fullfile(ModeFolder,'data_stressYZ.txt');   
    stressYZ = dlmread(DataFile);
    
    DataFile = fullfile(ModeFolder,'data_stressXZ.txt');   
    stressXZ = dlmread(DataFile);

end




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



    
    
    