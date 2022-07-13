
function [strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = get_ansys_strain(obj,mode_number)

    % Loads from the .txt files stored in output_files, for a given mode, the
    % element-by-element results for the strains, as calculated by Ansys.
    % The Cartesian reference frame in which they are expressed is dependent on the element type, 
    % and may vary from one element to another (e.g. from one element belonging to the substrate
    % to one belonging to the coating), 
    % however it is guaranteed to coincide
    % with the frame for the stresses on the SAME element,
    % so that the scalar energy quantities can be readily evaluated.
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    % mode_number - integer. 
    % strainXX, strainYY, etc. - arrays Nx1 double (N = number of elements).

    OutputFolder = obj.calculation_settings.apdl_output_folder;
    ModeFolder = fullfile(OutputFolder,num2str(mode_number));
    
    DataFile = fullfile(ModeFolder,'data_strainXX.txt');   
    strainXX = dlmread(DataFile);
    
    DataFile = fullfile(ModeFolder,'data_strainYY.txt');   
    strainYY = dlmread(DataFile); 
    
    DataFile = fullfile(ModeFolder,'data_strainZZ.txt');   
    strainZZ = dlmread(DataFile);  
    
    DataFile = fullfile(ModeFolder,'data_strainXY.txt');   
    strainXY = dlmread(DataFile);    
    
    DataFile = fullfile(ModeFolder,'data_strainYZ.txt');   
    strainYZ = dlmread(DataFile); 
    
    DataFile = fullfile(ModeFolder,'data_strainXZ.txt');   
    strainXZ = dlmread(DataFile);

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


