
function obj = load_elastic_energies(obj)

    % Loads into calculation_results a vector containing the results for elastic_energy,
    % for all modes included in mode_list, in the same order.
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).

    mode_list = obj.calculation_settings.mode_list;
    elastic_energy_modes = [];
    
    for mode_number = mode_list
        
        elastic_energy_modes(end+1) = sum(obj.get_element_elastic_energy(mode_number));

    end

    obj.calculation_results.elastic_energy_modes = elastic_energy_modes;

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



