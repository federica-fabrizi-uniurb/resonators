
function [positionsX, positionsY, positionsZ] = get_ansys_positions(obj,mode_number)

    % Loads from the .txt files stored in output_files, for a given mode, the
    % element-by-element results for the element position, in meters, expressed
    % in the Ansys reference frame (oriented like the Geometry Reference Frame used by this program).
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    % mode_number - integer. 
    % positionsX, positionsY, positionsZ - arrays Nx1 double (N = number of elements).

    OutputFolder = obj.calculation_settings.apdl_output_folder;
    ModeFolder = fullfile(OutputFolder,num2str(mode_number));
    
    DataFile = fullfile(ModeFolder,'data_positionsX.txt');   
    positionsX = dlmread(DataFile);
    
    DataFile = fullfile(ModeFolder,'data_positionsY.txt');   
    positionsY = dlmread(DataFile); 
    
    DataFile = fullfile(ModeFolder,'data_positionsZ.txt');   
    positionsZ = dlmread(DataFile);  
    
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



