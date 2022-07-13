
classdef Resonator 
    
    % Superclass for all resonator types (type Disc, DoublyCoatedDisc,
    % CurvedDisc, CantileverFiber are the subclasses implemented so far).
    
    properties
        
        substrate
        mesh_settings      
        calculation_settings
        calculation_results
        percentage_uncertainty_meas
        
    end
    
    properties (Hidden = true)
        
        path_ansys_exe = 'C:\Program Files\ANSYS Inc\v202\ansys\bin\winx64\ANSYS202.exe';
        path_workbench_exe = 'C:\Program Files\ANSYS Inc\v202\Framework\bin\Win64\RunWB2.exe';
        
    end
    
    methods
        
        function obj = Resonator()
            
            % constructor function, automatically called at instantiation
            
            % classes within classes: do not initialise properties in the class definition but do it here in the constructor
            % else the values returned by these expressions are part of the class definition and are constant for all instances of the class
            
            %obj.coating.diameter = obj.substrate.diameter;
            
            obj.substrate = Body();
            obj.percentage_uncertainty_meas = 1;
            obj.mesh_settings = MeshSettings();        
            obj.calculation_settings = CalculationSettings(class(obj));
            
        end
        
        function obj = set_folder(obj,object_folder)
            
            obj = local_set_folder(obj,object_folder);
            
        end
        
        function volume = get_element_volume(obj,mode_number)
            
            volume = get_ansys_volume(obj,mode_number);
            
        end
        
        function volumechange = get_element_volumechange(obj,mode_number)
            
            volumechange = get_ansys_volumechange(obj,mode_number);
            
        end
        
        function [strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = get_element_strain(obj,mode_number)
            
            [strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = get_ansys_strain(obj,mode_number);
            
        end

        function [stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = get_element_stress(obj,mode_number)
            
            [stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = get_ansys_stress(obj,mode_number);
            
        end
        
        function elastic_energy = get_element_elastic_energy(obj,mode_number)
            
            elastic_energy = get_ansys_sene(obj,mode_number);
            
        end
        
        function [elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = get_element_energy_approx(obj,mode_number)
            
            [elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = get_ansys_energy_approx(obj,mode_number);
            
        end 
        
        function save(obj)
            
            resonator = obj;
            save(fullfile(obj.calculation_settings.matlab_variable_folder,'resonator.mat'),'resonator');
            % to be re-loaded with:
            % myvar = load('path\to\object_folder\matlab_variable\resonator.mat').resonator;
            
        end
        
        function obj = reinitialise(obj)
            
            % delete object folder with all its data
            object_folder = obj.calculation_settings.object_folder;
            if exist(object_folder,'dir')
                   rmdir(object_folder,'s');
            end
            
            % delete matlab variable with all its properties and re-instantiate it
            constructor = class(obj);
            obj = feval(constructor);
            
            % re-create empty object folder (with subfolders)
            obj = obj.set_folder(object_folder);
            
        end    

    end
    
    methods (Access = protected)
        
        function obj = load_ansys_results(obj)            
            
            obj = load_frequencies(obj);
            obj = load_elastic_energies(obj);
            obj = load_energies_approx(obj);
            
        end
        
    end
    
    methods (Abstract)
        
        obj = set_defaults(obj)
        
        workbench_out = ansys_workbench_calculate(obj)
        
        apdl_out = ansys_apdl_calculate(obj) 
        
    end
    
end

% Local functions (only callable from within this same .m file)

function obj = local_set_folder(obj,object_folder)

    % create folder if not existing
    if ~exist(object_folder,'dir')
           mkdir(object_folder);
    end

    % if folder was given as a relative path, convert it to absolute  
    object_folder = what(object_folder).path;

    obj.calculation_settings.object_folder = object_folder;
    
    obj.calculation_settings.support_folder = fullfile(object_folder,'support_files');
    
    obj.calculation_settings.apdl_input_folder = fullfile(object_folder,'input_files');
    obj.calculation_settings.apdl_output_folder = fullfile(object_folder,'output_files');
    
    obj.calculation_settings.ansys_folder = fullfile(object_folder,'ansys_files');
    
    obj.calculation_settings.matlab_variable_folder = fullfile(object_folder,'matlab_variable');
    
    % create subfolders if not existing
    if ~exist(obj.calculation_settings.support_folder,'dir')
           mkdir(obj.calculation_settings.support_folder);
    end
    if ~exist(obj.calculation_settings.apdl_input_folder,'dir')
           mkdir(obj.calculation_settings.apdl_input_folder);
    end
    if ~exist(obj.calculation_settings.apdl_output_folder,'dir')
           mkdir(obj.calculation_settings.apdl_output_folder);
    end
    if ~exist(obj.calculation_settings.ansys_folder,'dir')
           mkdir(obj.calculation_settings.ansys_folder);
    end
    if ~exist(obj.calculation_settings.matlab_variable_folder,'dir')
       mkdir(obj.calculation_settings.matlab_variable_folder);
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



