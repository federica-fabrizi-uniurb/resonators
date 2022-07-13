
classdef CurvedDisc < Resonator
    
    % Subclass of Resonator, for type CurvedDisc;
    % contains methods to execute FEA with Ansys. 
    
    properties
        
        curvature_radius    % curvature_radius at half thickness:
                            % curvature_radius_external = curvature_radius + thickness/2
                            % curvature_radius_internal = curvature_radius - thickness/2
        
    end
    
    methods
        
        function obj = CurvedDisc()
            
            obj.calculation_results = CalculationResults_CurvedDisc();
            
        end
        
        function obj = set_defaults(obj)
            
            obj = local_set_defaults(obj);
            
        end
        
        function obj = set_curvature(obj,curvature_radius)
            
            obj = local_set_curvature(obj,curvature_radius);
            
        end
        
        function workbench_out = ansys_workbench_calculate(obj)
            
            workbench_prepare_input_file_CurvedDisc(obj);
            workbench_out = workbench_calculate(obj);
            
        end
        
        function obj = ansys_apdl_calculate(obj)    
            
            apdl_prepare_input_file_CurvedDisc(obj);
            apdl_out = apdl_calculate(obj);
            obj = load_ansys_results(obj);
            
        end
        
    end
    
    methods (Access = protected)
        
        function obj = load_ansys_results(obj)            
            
            obj = load_ansys_results@Resonator(obj);
            obj = local_load_D_TED(obj);
            
        end
        
    end
    
end

% Local functions (only callable from within this same .m file)

function obj = local_set_defaults(obj)

    obj.substrate.diameter = 0.0254;
    obj.substrate.thickness = .3e-3;

    obj = obj.set_curvature(1);  % Curvature radius (supposed to be large)

    obj.substrate.material = Material('silica');
    obj.substrate.material = obj.substrate.material.set_temperature(300);

end

function obj = local_set_curvature(obj,curvature_radius)
    
    obj.curvature_radius = curvature_radius;
    
    thickness = obj.substrate.thickness;
    diameter = obj.substrate.diameter;

    ext_radius = curvature_radius + thickness/2;
    int_radius = curvature_radius - thickness/2;

    Omega = 3*pi*(diameter/2)^2*(ext_radius^2 + ext_radius*int_radius + int_radius^2)^-1;
    cos_theta = 1 - Omega/(2*pi);
    theta = acos(cos_theta);
    
    volume_curved = Omega/3*(ext_radius^3-int_radius^3);
    volume_flat = pi*(diameter/2)^2*thickness;

    obj.calculation_results.volume_curved = volume_curved;
    obj.calculation_results.volume_flat = volume_flat;
    obj.calculation_results.curvature_solid_angle = Omega;
    obj.calculation_results.curvature_planar_angle = 2*theta;  % full span (2*theta)

end

function obj = local_load_D_TED(obj)

    obj.calculation_results.D_TED_modes = obj.calculation_results.dilatation_energy_approx_modes./obj.calculation_results.elastic_energy_approx_modes;

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



