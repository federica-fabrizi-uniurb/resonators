
classdef DoublyCoatedDisc < Resonator
    
    % Subclass of Resonator, for type DoublyCoatedDisc;
    % contains methods to execute FEA with Ansys and 
    % to calculate the thermoelastic dissipation analytically.   
    
    properties
        
        phi_coating
        coating
        %%%heat_flow_direction = [0 0 1]   % in sysref attached to the geometry (axis 3 = axis of highest symmetry)
        
    end
    
    methods
        
        function obj = DoublyCoatedDisc()
            
            obj.coating = Body();
            obj.phi_coating = 1e-4;
            obj.calculation_results = CalculationResults_DoublyCoatedDisc();
            
        end
        
        function obj = set_defaults(obj)
            
            obj = local_set_defaults(obj);
            
        end
        
        function workbench_out = ansys_workbench_calculate(obj)
            
            workbench_prepare_input_file_DoublyCoatedDisc(obj);
            workbench_out = workbench_calculate(obj);
            
        end
        
        function obj = ansys_apdl_calculate(obj)    
            
            apdl_prepare_input_file_DoublyCoatedDisc(obj);
            apdl_out = apdl_calculate(obj);
            obj = load_ansys_results(obj);
            
        end
        
        function obj = calculate_TED_vengallatore(obj,n_veng)

            if nargin < 2
                n_veng = 1;
            end
            obj = local_calculate_TED_vengallatore(obj,n_veng);

        end  
        
        function obj = calculate_dissipation_TED_plus_coating(obj)
            
            obj = local_calculate_dissipation_TED_plus_coating(obj);
            
        end
        
    end
    
    methods (Access = protected)
        
        function obj = load_ansys_results(obj)            
            
            obj = load_ansys_results@Resonator(obj);
            obj = load_element_numbers(obj);
            obj = load_Dc(obj);
            obj = local_load_D_TED(obj);
            
        end 
        
    end
    
end

% Local functions (only callable from within this same .m file)

function obj = local_set_defaults(obj)

    obj.substrate.diameter = 0.0254*2;
    obj.substrate.thickness = .2e-3;

    obj.coating.thickness = 1e-6;

    obj.substrate.material = Material('silicon');
    obj.substrate.material.flag_anisotropic = 'single_crystal';
    obj.substrate.material.orientation_axis_1 = [0 0 1];	

    obj.coating.material = Material('silica');

    obj.substrate.material = obj.substrate.material.set_temperature(300);
    obj.coating.material = obj.coating.material.set_temperature(300);


end

function obj = local_load_D_TED(obj)

    obj.calculation_results.D_TED_modes = obj.calculation_results.dilatation_energy_approx_substrate_modes./obj.calculation_results.elastic_energy_approx_substrate_modes;

end

function obj = local_calculate_TED_vengallatore(obj,n_veng)
    
    temperature = obj.substrate.material.temperature;
    
    % Create material properties structure to pass to veng_c
    params = [];
    
    params.Cs = obj.substrate.material.specific_heat_capacity*obj.substrate.material.density;
    % if material is amorphous (isotropic)
    if strcmp(obj.substrate.material.flag_anisotropic,'amorph') 
        % thermal conductivity
        if ~isempty(obj.substrate.material.thermal_conductivity_amorph)
            params.ks = obj.substrate.material.thermal_conductivity_amorph;
        else
            params.ks = obj.substrate.material.thermal_conductivity_poly;
        end
        % thermal linear expansion
        if ~isempty(obj.substrate.material.thermal_linear_expansion_amorph)
            params.alphas = obj.substrate.material.thermal_linear_expansion_amorph;
        else
            params.alphas = obj.substrate.material.thermal_linear_expansion_poly;
        end    
        % young
        if ~isempty(obj.substrate.material.young_amorph)
            params.Es = obj.substrate.material.young_amorph;
        else
            params.Es = obj.substrate.material.young_poly;
        end
        % bulk
        if ~isempty(obj.substrate.material.bulk_amorph)
            params.BMs = obj.substrate.material.bulk_amorph;
        else
            if ~isempty(obj.substrate.material.young_amorph) && ~isempty(obj.substrate.material.poisson_amorph)
                params.BMs = obj.substrate.material.bulk_amorph_from_young_amorph_and_poisson_amorph();
            else
                if ~isempty(obj.substrate.material.bulk_poly)
                    params.BMs = obj.substrate.material.bulk_poly;
                else
                    params.BMs = obj.substrate.material.bulk_poly_from_young_poly_and_poisson_poly();
                end
            end
        end 
    else
        % if material is polycrystalline (isotropic)    
        if strcmp(obj.substrate.material.flag_anisotropic,'poly') 
            % thermal conductivity
            if ~isempty(obj.substrate.material.thermal_conductivity_poly)
                params.ks = obj.substrate.material.thermal_conductivity_poly;
            else
                params.ks = 1/3*trace(obj.substrate.material.thermal_conductivity_matrix);
            end
            % thermal linear expansion
            if ~isempty(obj.substrate.material.thermal_linear_expansion_poly)
                params.alphas = obj.substrate.material.thermal_linear_expansion_poly;
            else
                params.alphas = obj.substrate.material.alpha_poly_from_alpha_matrix_and_stiffness_matrix();
            end             
            % young
            if ~isempty(obj.substrate.material.young_poly)
                params.Es = obj.substrate.material.young_poly;
            else
                params.Es = obj.substrate.material.young_poly_from_stiffness_matrix();
            end
            % bulk
            if ~isempty(obj.substrate.material.bulk_poly)
                params.BMs = obj.substrate.material.bulk_poly;
            else
                if ~isempty(obj.substrate.material.young_poly) && ~isempty(obj.substrate.material.poisson_poly)
                    params.BMs = obj.substrate.material.bulk_poly_from_young_poly_and_poisson_poly();
                else
                    params.BMs = obj.substrate.material.bulk_poly_from_stiffness_matrix();
                end
            end
        else
            % if material is single crystal (anisotropic)  
            if strcmp(obj.substrate.material.flag_anisotropic,'single_crystal') 
                %%% hf_direction = obj.heat_flow_direction; in MRF
                hf_direction = obj.substrate.material.orientation_axis_1;
                [hf_p1_direction, hf_p2_direction] = find_orthogonal_vectors(hf_direction);
                % thermal conductivity
                if ~isempty(obj.substrate.material.thermal_conductivity_matrix)
                    params.ks = obj.substrate.material.k_directional_from_k_matrix(hf_direction);
                else
                    params.ks = obj.substrate.material.thermal_conductivity_directional.(strrep(['direction_' num2str(hf_direction(1)) num2str(hf_direction(2)) num2str(hf_direction(3))],'-','m'));
                end
                % thermal linear expansion
                if ~isempty(obj.substrate.material.thermal_linear_expansion_matrix)
                    params.alphas = mean([obj.substrate.material.alpha_directional_from_alpha_matrix(hf_p1_direction),obj.substrate.material.alpha_directional_from_alpha_matrix(hf_p2_direction)]);
                else
                    alpha_perp = [];
                    alpha_dirs =  fieldnames(obj.substrate.material.thermal_linear_expansion_directional);
                    for jj = 1:length(alpha_dirs)
                        alpha_dir_string = erase(alpha_dirs{jj},'direction_');
                        alpha_dir_string = replace(alpha_dir_string,'m','-');
                        alpha_dir = str2double(regexp(alpha_dir_string,'(-|)(\d)','match')');
                        if abs(sum(alpha_dir.*hf_direction)) < 0.01 
                            alpha_perp(end+1) = obj.substrate.material.thermal_linear_expansion_directional.(alpha_dirs{jj});
                        end
                    end
                    params.alphas = mean(alpha_perp);                 
                end     
                % young
                if ~isempty(obj.substrate.material.stiffness_matrix)
                    params.Es = mean([obj.substrate.material.young_directional_from_stiffness_matrix(hf_p1_direction),obj.substrate.material.young_directional_from_stiffness_matrix(hf_p2_direction)]);
                else
                    E_perp = [];
                    E_dirs =  fieldnames(obj.substrate.material.young_directional);
                    for jj = 1:length(E_dirs)
                        E_dir_string = erase(E_dirs{jj},'direction_');
                        E_dir_string = replace(E_dir_string,'m','-');
                        E_dir = str2double(regexp(E_dir_string,'(-|)(\d)','match')');
                        if abs(sum(E_dir.*hf_direction)) < 0.01 
                            E_perp(end+1) = obj.substrate.material.young_directional.(E_dirs{jj});
                        end
                    end
                    params.Es = mean(E_perp);
                end
                % bulk
                if ~isempty(obj.substrate.material.stiffness_matrix)
                    params.BMs = obj.substrate.material.bulk_single_crystal_from_stiffness_matrix();
                else
                    young_mean = mean(cell2mat(struct2cell(obj.substrate.material.young_directional)));
                    poisson_mean = mean(cell2mat(struct2cell(obj.substrate.material.poisson_directional)));
                    params.BMs = young_mean/(3-6*poisson_mean);
                end                
            end
        end
    end    
    
    params.Cc = obj.coating.material.specific_heat_capacity*obj.coating.material.density;
    % if material is amorphous (isotropic)
    if strcmp(obj.coating.material.flag_anisotropic,'amorph') 
        % thermal conductivity
        if ~isempty(obj.coating.material.thermal_conductivity_amorph)
            params.kc = obj.coating.material.thermal_conductivity_amorph;
        else
            params.kc = obj.coating.material.thermal_conductivity_poly;
        end
        % thermal linear expansion
        if ~isempty(obj.coating.material.thermal_linear_expansion_amorph)
            params.alphac = obj.coating.material.thermal_linear_expansion_amorph;
        else
            params.alphac = obj.coating.material.thermal_linear_expansion_poly;
        end    
        % young
        if ~isempty(obj.coating.material.young_amorph)
            params.Ec = obj.coating.material.young_amorph;
        else
            params.Ec = obj.coating.material.young_poly;
        end
    else
        % if material is polycrystalline (isotropic)    
        if strcmp(obj.coating.material.flag_anisotropic,'poly') 
            % thermal conductivity
            if ~isempty(obj.coating.material.thermal_conductivity_poly)
                params.kc = obj.coating.material.thermal_conductivity_poly;
            else
                params.kc = 1/3*trace(obj.coating.material.thermal_conductivity_matrix);
            end
            % thermal linear expansion
            if ~isempty(obj.coating.material.thermal_linear_expansion_poly)
                params.alphac = obj.coating.material.thermal_linear_expansion_poly;
            else
                params.alphac = obj.coating.material.alpha_poly_from_alpha_matrix_and_stiffness_matrix();
            end             
            % young
            if ~isempty(obj.coating.material.young_poly)
                params.Ec = obj.coating.material.young_poly;
            else
                params.Ec = obj.coating.material.young_poly_from_stiffness_matrix();
            end
        else
            % if material is single crystal (anisotropic)  
            if strcmp(obj.coating.material.flag_anisotropic,'single_crystal') 
                %%% hf_direction = obj.heat_flow_direction;
                hf_direction = obj.coating.material.orientation_axis_1;
                [hf_p1_direction, hf_p2_direction] = find_orthogonal_vectors(hf_direction);
                % thermal conductivity
                if ~isempty(obj.coating.material.thermal_conductivity_matrix)
                    params.kc = obj.coating.material.k_directional_from_k_matrix(hf_direction);
                else
                    params.kc = obj.coating.material.thermal_conductivity_directional.(strrep(['direction_' num2str(hf_direction(1)) num2str(hf_direction(2)) num2str(hf_direction(3))],'-','m'));
                end
                % thermal linear expansion
                if ~isempty(obj.coating.material.thermal_linear_expansion_matrix)
                    params.alphac = mean([obj.coating.material.alpha_directional_from_alpha_matrix(hf_p1_direction),obj.coating.material.alpha_directional_from_alpha_matrix(hf_p2_direction)]);
                else
                    alpha_perp = [];
                    alpha_dirs =  fieldnames(obj.coating.material.thermal_linear_expansion_directional);
                    for jj = 1:length(alpha_dirs)
                        alpha_dir_string = erase(alpha_dirs{jj},'direction_');
                        alpha_dir_string = replace(alpha_dir_string,'m','-');
                        alpha_dir = str2double(regexp(alpha_dir_string,'(-|)(\d)','match')');
                        if abs(sum(alpha_dir.*hf_direction)) < 0.01 
                            alpha_perp(end+1) = obj.coating.material.thermal_linear_expansion_directional.(alpha_dirs{jj});
                        end
                    end
                    params.alphac = mean(alpha_perp);
                end     
                % young
                if ~isempty(obj.coating.material.stiffness_matrix)
                    params.Ec = mean([obj.coating.material.young_directional_from_stiffness_matrix(hf_p1_direction),obj.coating.material.young_directional_from_stiffness_matrix(hf_p2_direction)]);
                else
                    E_perp = [];
                    E_dirs =  fieldnames(obj.coating.material.young_directional);
                    for jj = 1:length(E_dirs)
                        E_dir_string = erase(E_dirs{jj},'direction_');
                        E_dir_string = replace(E_dir_string,'m','-');
                        E_dir = str2double(regexp(E_dir_string,'(-|)(\d)','match')');
                        if abs(sum(E_dir.*hf_direction)) < 0.01 
                            E_perp(end+1) = obj.coating.material.young_directional.(E_dirs{jj});
                        end
                    end
                    params.Ec = mean(E_perp);
                end               
            end
        end
    end    
    
    params.SUB = obj.substrate.material.keyword;
    params.COAT = obj.coating.material.keyword;

    % calculate
    a = obj.substrate.thickness/2;
    b = obj.substrate.thickness/2 + obj.coating.thickness; 

    %[ff,phiTED] = SVthel2(0,obj.coating.thickness,a,0,0,n_veng,params,temperature); % WITH COATING
    [ff,phiTED] = veng_c(a,b,n_veng,params,temperature); % WITH COATING

    phi_TED_veng_modes = [];
    for frequency = obj.calculation_results.frequency_modes
        phi_TED_veng_modes(end+1) = interp1(ff,phiTED,frequency);
    end

    % correction for D_TED
    D_TED_Vengallatore = params.Es/(9*params.BMs);
    phi_TED_without_D_TED_modes = phi_TED_veng_modes/D_TED_Vengallatore;
    phi_TED_modes = phi_TED_without_D_TED_modes.*obj.calculation_results.D_TED_modes;
    % note that we doe NOT include 2*pi factor in Zener's modulus
    
    ff_phi_TED_without_D_TED = phiTED/D_TED_Vengallatore;
    
    % load results
    obj.calculation_results.phi_TED_modes = phi_TED_modes;
    obj.calculation_results.phi_TED_without_D_TED_modes = phi_TED_without_D_TED_modes;
    obj.calculation_results.ff = ff;
    obj.calculation_results.ff_phi_TED_without_D_TED = ff_phi_TED_without_D_TED;
    
end

function  obj = local_calculate_dissipation_TED_plus_coating(obj)

    phi_MEAS_modes = obj.calculation_results.phi_TED_modes + obj.calculation_results.D_c_modes*obj.phi_coating;
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


