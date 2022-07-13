
function obj = load_Dc(obj)

    % Loads into calculation_results 12 vectors containing the results for the 4 quantities
    % elastic_energy, elastic_energy_approx, dilatation_energy_approx, and shear_energy_approx, 
    % each separated into its 3 cntributions from substrate, coating1 and coating2. 
    % Finally, loads into calculation_results a 13th vector, containing the
    % dilution factors D_c for the coatings, calculated from the results above.
    % Each vector stores the values for all modes included in mode_list, in the same order as they appear in the list.
    % This function is only usable for a resonator of type DoublyCoatedDisc
    % (as opposed to functions load_elastic_energies and load_energies_approx, which load the total results without 
    % distinguishing between the elements, and therefore are used for a resonator of any type).
    %
    % obj - instance of class DoublyCoatedDisc.
    
    mode_list = obj.calculation_settings.mode_list;
    element_number_substrate = obj.calculation_results.element_number_substrate;
    element_number_coating1 = obj.calculation_results.element_number_coating1;
    element_number_coating2 = obj.calculation_results.element_number_coating2;
    
    elastic_energy_substrate_modes = [];
    elastic_energy_coating1_modes = [];
    elastic_energy_coating2_modes = [];

    elastic_energy_approx_substrate_modes = [];
    elastic_energy_approx_coating1_modes = [];
    elastic_energy_approx_coating2_modes = [];     
    
    dilatation_energy_approx_substrate_modes = [];
    dilatation_energy_approx_coating1_modes = [];
    dilatation_energy_approx_coating2_modes = [];
    
    shear_energy_approx_substrate_modes = [];
    shear_energy_approx_coating1_modes = [];
    shear_energy_approx_coating2_modes = [];    
    
    for mode_number = mode_list
        
        element_elastic_energy = obj.get_element_elastic_energy(mode_number);
        
        elastic_energy_substrate_modes(end+1) = sum(element_elastic_energy(1:element_number_substrate));
        elastic_energy_coating1_modes(end+1) = sum(element_elastic_energy(element_number_substrate+1:element_number_substrate+element_number_coating1));
        elastic_energy_coating2_modes(end+1) = sum(element_elastic_energy(element_number_substrate+element_number_coating1+1:element_number_substrate+element_number_coating1+element_number_coating2));
    
        [element_elastic_energy_approx,element_dilatation_energy_approx,element_shear_energy_approx] = obj.get_element_energy_approx(mode_number);
        
        elastic_energy_approx_substrate_modes(end+1) = sum(element_elastic_energy_approx(1:element_number_substrate));
        elastic_energy_approx_coating1_modes(end+1) = sum(element_elastic_energy_approx(element_number_substrate+1:element_number_substrate+element_number_coating1));
        elastic_energy_approx_coating2_modes(end+1) = sum(element_elastic_energy_approx(element_number_substrate+element_number_coating1+1:element_number_substrate+element_number_coating1+element_number_coating2));
    
        dilatation_energy_approx_substrate_modes(end+1) = sum(element_dilatation_energy_approx(1:element_number_substrate));
        dilatation_energy_approx_coating1_modes(end+1) = sum(element_dilatation_energy_approx(element_number_substrate+1:element_number_substrate+element_number_coating1));
        dilatation_energy_approx_coating2_modes(end+1) = sum(element_dilatation_energy_approx(element_number_substrate+element_number_coating1+1:element_number_substrate+element_number_coating1+element_number_coating2));
  
        shear_energy_approx_substrate_modes(end+1) = sum(element_shear_energy_approx(1:element_number_substrate));
        shear_energy_approx_coating1_modes(end+1) = sum(element_shear_energy_approx(element_number_substrate+1:element_number_substrate+element_number_coating1));
        shear_energy_approx_coating2_modes(end+1) = sum(element_shear_energy_approx(element_number_substrate+element_number_coating1+1:element_number_substrate+element_number_coating1+element_number_coating2));
        
    end

    obj.calculation_results.elastic_energy_substrate_modes = elastic_energy_substrate_modes;
    obj.calculation_results.elastic_energy_coating1_modes = elastic_energy_coating1_modes;
    obj.calculation_results.elastic_energy_coating2_modes = elastic_energy_coating2_modes;
    
    obj.calculation_results.elastic_energy_approx_substrate_modes = elastic_energy_approx_substrate_modes;
    obj.calculation_results.elastic_energy_approx_coating1_modes = elastic_energy_approx_coating1_modes;
    obj.calculation_results.elastic_energy_approx_coating2_modes = elastic_energy_approx_coating2_modes;    
    
    obj.calculation_results.dilatation_energy_approx_substrate_modes = dilatation_energy_approx_substrate_modes;
    obj.calculation_results.dilatation_energy_approx_coating1_modes = dilatation_energy_approx_coating1_modes;
    obj.calculation_results.dilatation_energy_approx_coating2_modes = dilatation_energy_approx_coating2_modes;      
    
    obj.calculation_results.shear_energy_approx_substrate_modes = shear_energy_approx_substrate_modes;
    obj.calculation_results.shear_energy_approx_coating1_modes = shear_energy_approx_coating1_modes;
    obj.calculation_results.shear_energy_approx_coating2_modes = shear_energy_approx_coating2_modes;     
    
    obj.calculation_results.D_c_modes = (elastic_energy_coating1_modes + elastic_energy_coating2_modes)./(elastic_energy_substrate_modes + elastic_energy_coating1_modes + elastic_energy_coating2_modes);

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

