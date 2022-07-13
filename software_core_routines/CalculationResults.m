
classdef CalculationResults 
    
    % Class containing the calculation results common to every resonator type;
    % this class features subclasses adding further results, specific only to some resonator types.
    
    properties
       
        frequency_modes
        elastic_energy_modes
        elastic_energy_approx_modes
        dilatation_energy_approx_modes
        shear_energy_approx_modes
        D_TED_modes
        ff
        ff_phi_TED_without_D_TED
        phi_TED_modes                 % with_D_TED
        phi_TED_without_D_TED_modes
        phi_MEAS_modes
        Delta_phi_MEAS_modes        
        
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


