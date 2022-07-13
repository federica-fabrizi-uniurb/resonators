
classdef CalculationSettings 
    
    % Class containing the settings for the Ansys calculation, requested by the user.    
    
    properties
       
        mode_list
        min_frequency
        max_frequency
        max_number_of_modes
        object_folder
        
    end
    
    properties %(Hidden=true)
        
        templates_folder 
        workbench_input_file_journal_template
        workbench_input_file_geometry_template
        workbench_input_file_mechanical_template
        
        support_folder
        workbench_input_file_journal_edited
        workbench_input_file_geometry_edited
        workbench_input_file_mechanical_edited
        
        apdl_input_folder 
        apdl_input_file_template
        apdl_input_file_edited
        
        apdl_output_folder
        
        ansys_folder
        
        matlab_variable_folder

    end
    
    methods
        
        function obj = CalculationSettings(object_type)
            
            % constructor function, automatically called at instantiation
            
            obj.mode_list = [1:10];
            obj.min_frequency = 100;
            obj.max_frequency = 100000;
            obj.max_number_of_modes = 30;
        
     
            obj.templates_folder = fullfile(fileparts(mfilename('fullpath')),['templates_',object_type]);
            
            obj.workbench_input_file_journal_template = 'journal.wbjn';
            obj.workbench_input_file_geometry_template = 'geom.js';
            obj.workbench_input_file_mechanical_template = 'mech.py';
            
            obj.workbench_input_file_journal_edited = 'journal_edited.wbjn';
            obj.workbench_input_file_geometry_edited = 'geom_edited.js';
            obj.workbench_input_file_mechanical_edited = 'mech_edited.py';
            
            obj.apdl_input_file_template = 'input_file.dat';
            obj.apdl_input_file_edited = 'input_file_edited.dat';
            
        end
        
    end
    
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



