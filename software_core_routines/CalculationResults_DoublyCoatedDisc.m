
classdef CalculationResults_DoublyCoatedDisc < CalculationResults
    
    % Subclass of CalculationResults, adding the results specific only to a resonator of type DoublyCoatedDisc.
    
    properties
       
        element_number_substrate
        element_number_coating1
        element_number_coating2
        
        elastic_energy_substrate_modes
        elastic_energy_coating1_modes
        elastic_energy_coating2_modes
        
        elastic_energy_approx_substrate_modes
        elastic_energy_approx_coating1_modes
        elastic_energy_approx_coating2_modes  
        
        dilatation_energy_approx_substrate_modes
        dilatation_energy_approx_coating1_modes
        dilatation_energy_approx_coating2_modes         
        
        shear_energy_approx_substrate_modes
        shear_energy_approx_coating1_modes
        shear_energy_approx_coating2_modes       
        
        D_c_modes
        
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


