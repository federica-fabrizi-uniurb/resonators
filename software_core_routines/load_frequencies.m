
function obj = load_frequencies(obj)

    % Loads into calculation_results a vector containing the modal frequencies calculated by Ansys,
    % for all modes included in mode_list, in the same order. 
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    
    any_mode_number = obj.calculation_settings.mode_list(1);
    OutputFolder = obj.calculation_settings.apdl_output_folder;
    ModeFolder = fullfile(OutputFolder,num2str(any_mode_number));

    DataFile = fullfile(ModeFolder,'data_frequencies.txt');
    
    data = dlmread(DataFile);
    
    frequency_all_modes = data(:,1);
    
    % keep only the modes in mode_list
    
    mode_list = obj.calculation_settings.mode_list;
    frequency_modes = [];
    
    for mode_number = mode_list
        
        frequency_modes(end+1) = frequency_all_modes(mode_number);
        
    end
    
    obj.calculation_results.frequency_modes = frequency_modes;
    
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



    
    