
classdef CantileverFiber < Resonator
    
    % Subclass of Resonator, for type CantileverFiber;
    % contains methods to execute FEA with Ansys and 
    % to calculate the thermoelastic dissipation analytically.
    
    properties
        
        %%%heat_flow_direction = [1 0 0]   % in GRF, reference system attached to the geometry (perpendicular to axis Z = axis of highest symmetry)

    end
    
    methods
        
        function obj = CantileverFiber()
            
            obj.calculation_results = CalculationResults_Disc();
            obj.mesh_settings.divisions = 0;
            
        end
        
        function obj = set_defaults(obj)
            
            obj = local_set_defaults(obj);
            
        end
        
        function workbench_out = ansys_workbench_calculate(obj)
            
            workbench_prepare_input_file_Disc(obj);
            workbench_out = workbench_calculate(obj);
            
        end
        
        function obj = ansys_apdl_calculate(obj)    
            
            apdl_prepare_input_file_Disc(obj);
            apdl_out = apdl_calculate(obj);
            obj = load_ansys_results(obj);
            
        end
        
        function obj = calculate_TED_zener(obj)

            obj = local_calculate_TED_zener(obj);

        end  
        
        function obj = calculate_dissipation_TED(obj)
            
            obj = local_calculate_dissipation_TED(obj);
            
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

    obj.substrate.diameter = 2e-3;     % [m] 2 mm
    obj.substrate.thickness = 335e-3;  % [m] 335 mm

    obj.substrate.material = Material('silica');
    obj.substrate.material = obj.substrate.material.set_temperature(300);

end

function obj = local_load_D_TED(obj)

    obj.calculation_results.D_TED_modes = obj.calculation_results.dilatation_energy_approx_modes./obj.calculation_results.elastic_energy_approx_modes;

end

function obj = local_calculate_TED_zener(obj)

    % Calculation done according to 
    % 'Study of the Thermoelastic Effect in a Cylinder'
    % F. Piergiovanni, F. Martelli
    
    % get relevant properties
    
    temperature = obj.substrate.material.temperature;
    
    density = obj.substrate.material.density;
    specific_heat_capacity = obj.substrate.material.specific_heat_capacity;
    
    % if material is amorphous (isotropic)
    if strcmp(obj.substrate.material.flag_anisotropic,'amorph') 
        % thermal conductivity
        if ~isempty(obj.substrate.material.thermal_conductivity_amorph)
            thermal_conductivity = obj.substrate.material.thermal_conductivity_amorph;
        else
            thermal_conductivity = obj.substrate.material.thermal_conductivity_poly;
        end
        % thermal linear expansion
        if ~isempty(obj.substrate.material.thermal_linear_expansion_amorph)
            thermal_linear_expansion = obj.substrate.material.thermal_linear_expansion_amorph;
        else
            thermal_linear_expansion = obj.substrate.material.thermal_linear_expansion_poly;
        end    
        % young
        if ~isempty(obj.substrate.material.young_amorph)
            young = obj.substrate.material.young_amorph;
        else
            young = obj.substrate.material.young_poly;
        end
        % bulk
        if ~isempty(obj.substrate.material.bulk_amorph)
            bulk = obj.substrate.material.bulk_amorph;
        else
            if ~isempty(obj.substrate.material.young_amorph) && ~isempty(obj.substrate.material.poisson_amorph)
                bulk = obj.substrate.material.bulk_amorph_from_young_amorph_and_poisson_amorph();
            else
                if ~isempty(obj.substrate.material.bulk_poly)
                    bulk = obj.substrate.material.bulk_poly;
                else
                    bulk = obj.substrate.material.bulk_poly_from_young_poly_and_poisson_poly();
                end
            end
        end 
    else
        % if material is polycrystalline (isotropic)    
        if strcmp(obj.substrate.material.flag_anisotropic,'poly') 
            % thermal conductivity
            if ~isempty(obj.substrate.material.thermal_conductivity_poly)
                thermal_conductivity = obj.substrate.material.thermal_conductivity_poly;
            else
                thermal_conductivity = 1/3*trace(obj.substrate.material.thermal_conductivity_matrix);
            end
            % thermal linear expansion
            if ~isempty(obj.substrate.material.thermal_linear_expansion_poly)
                thermal_linear_expansion = obj.substrate.material.thermal_linear_expansion_poly;
            else
                thermal_linear_expansion = obj.substrate.material.alpha_poly_from_alpha_matrix();
            end             
            % young 
            if ~isempty(obj.substrate.material.young_poly)
                young = obj.substrate.material.young_poly;
            else
                young = obj.substrate.material.young_poly_from_stiffness_matrix();
            end
            % bulk
            if ~isempty(obj.substrate.material.bulk_poly)
                bulk = obj.substrate.material.bulk_poly;
            else
                if ~isempty(obj.substrate.material.young_poly) && ~isempty(obj.substrate.material.poisson_poly)
                    bulk = obj.substrate.material.bulk_poly_from_young_poly_and_poisson_poly();
                else
                    bulk = obj.substrate.material.bulk_poly_from_stiffness_matrix();
                end
            end
        else
            % if material is single crystal (anisotropic)  
            if strcmp(obj.substrate.material.flag_anisotropic,'single_crystal') 
                %%% hf_direction = obj.heat_flow_(direction; in MRF
                [hf_direction1, hf_direction2] = find_orthogonal_vectors(obj.substrate.material.orientation_axis_1);
                hf_p_direction = obj.substrate.material.orientation_axis_1;
                % thermal conductivity
                if ~isempty(obj.substrate.material.thermal_conductivity_matrix)
                    k1 = obj.substrate.material.k_directional_from_k_matrix(hf_direction1);
                    k2 = obj.substrate.material.k_directional_from_k_matrix(hf_direction2);
                    thermal_conductivity = mean([k1,k2]);
                else
                    k_perp = [];
                    k_dirs =  fieldnames(obj.substrate.material.thermal_conductivity_directional);
                    for jj = 1:length(k_dirs)
                        k_dir_string = erase(k_dirs{jj},'direction_');
                        k_dir_string = replace(k_dir_string,'m','-');
                        k_dir = str2double(regexp(k_dir_string,'(-|)(\d)','match')');
                        if abs(sum(k_dir.*hf_p_direction)) < 0.01 
                            k_perp(end+1) = obj.substrate.material.thermal_conductivity_directional.(k_dirs{jj});
                        end
                    end
                    thermal_conductivity = mean(k_perp);
                end
                % thermal linear expansion
                if ~isempty(obj.substrate.material.thermal_linear_expansion_matrix)
                    thermal_linear_expansion = obj.substrate.material.alpha_directional_from_alpha_matrix(hf_p_direction);
                else
                    alpha_dirs =  fieldnames(obj.substrate.material.thermal_linear_expansion_directional);
                    for jj = 1:length(alpha_dirs)
                        alpha_dir_string = erase(alpha_dirs{jj},'direction_');
                        alpha_dir_string = replace(alpha_dir_string,'m','-');
                        alpha_dir = str2double(regexp(alpha_dir_string,'(-|)(\d)','match')');
                        para_factors = alpha_dir./hf_p_direction;
                        if abs(para_factors(1) - para_factors(2)) < 0.01 && abs(para_factors(1) - para_factors(3)) < 0.01
                            alpha_para = obj.substrate.material.thermal_linear_expansion_directional.(alpha_dirs{jj});
                        end
                    end
                    thermal_linear_expansion = alpha_para;
                end     
                % young
                if ~isempty(obj.substrate.material.stiffness_matrix)
                    young = obj.substrate.material.young_directional_from_stiffness_matrix(hf_p_direction);
                else
                    E_dirs =  fieldnames(obj.substrate.material.young_directional);
                    for jj = 1:length(E_dirs)
                        E_dir_string = erase(E_dirs{jj},'direction_');
                        E_dir_string = replace(E_dir_string,'m','-');
                        E_dir = str2double(regexp(E_dir_string,'(-|)(\d)','match')');
                        para_factors = E_dir./hf_p_direction;
                        if abs(para_factors(1) - para_factors(2)) < 0.01 && abs(para_factors(1) - para_factors(3)) < 0.01
                            E_para = obj.substrate.material.young_directional.(E_dirs{jj});
                        end
                    end
                    young = E_para;                        
                end
                % bulk
                if ~isempty(obj.substrate.material.stiffness_matrix)
                    bulk = obj.substrate.material.bulk_single_crystal_from_stiffness_matrix();
                else
                    young_mean = mean(cell2mat(struct2cell(obj.substrate.material.young_directional)));
                    poisson_mean = mean(cell2mat(struct2cell(obj.substrate.material.poisson_directional)));
                    bulk = young_mean/(3-6*poisson_mean);
                end                
            end
        end
    end
    
    diameter = obj.substrate.diameter;

    % calculate
    
    ff = logspace(0,6,5000);
    
    zener_coeff = young*thermal_linear_expansion^2*temperature/(specific_heat_capacity*density);    
    omega_peak = (2*1.841)^2*thermal_conductivity/(specific_heat_capacity*density*diameter^2);
    
    phiTED = zener_coeff*2*pi*ff*omega_peak./((2*pi*ff).^2+omega_peak^2);

    phi_TED_zener_modes = [];
    for frequency = obj.calculation_results.frequency_modes
        phi_TED_zener_modes(end+1) = interp1(ff,phiTED,frequency);
    end
    % note that we do NOT include 2*pi factor in Zener modulus

    % correction for D_TED: 
    % in the case of a fiber,
    % D_TED calculated by Ansys should be the same for the various modes
    % and approx equal to young/(9*bulk)
    
    ff_phi_TED_without_D_TED = phiTED/young*(9*bulk);
    
    % load results
    obj.calculation_results.phi_TED_modes = phi_TED_zener_modes;
    obj.calculation_results.phi_TED_without_D_TED_modes = phi_TED_zener_modes./obj.calculation_results.D_TED_modes;
    obj.calculation_results.ff = ff;
    obj.calculation_results.ff_phi_TED_without_D_TED = ff_phi_TED_without_D_TED;
    
end

function  obj = local_calculate_dissipation_TED(obj)

    phi_MEAS_modes = obj.calculation_results.phi_TED_modes;
    obj.calculation_results.phi_MEAS_modes = phi_MEAS_modes;

    Delta_phi_MEAS_modes = phi_MEAS_modes/100*obj.percentage_uncertainty_meas;
    obj.calculation_results.Delta_phi_MEAS_modes = Delta_phi_MEAS_modes;

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

