
function obj = load_energies_approx(obj)

    % Loads into calculation_results three vectors containing the results for
    % elastic_energy_approx, dilatation_energy_approx, and shear_energy_approx, respectively,
    % for all modes included in mode_list, in the same order.
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).

    mode_list = obj.calculation_settings.mode_list;
    elastic_energy_approx_modes = [];
    dilatation_energy_approx_modes = [];
    shear_energy_approx_modes = [];
    
    for mode_number = mode_list
        
        [elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = obj.get_element_energy_approx(mode_number);
        elastic_energy_approx_modes(end+1) = sum(elastic_energy_approx);
        dilatation_energy_approx_modes(end+1) = sum(dilatation_energy_approx);
        shear_energy_approx_modes(end+1) = sum(shear_energy_approx);

    end

    obj.calculation_results.elastic_energy_approx_modes = elastic_energy_approx_modes;
    obj.calculation_results.dilatation_energy_approx_modes = dilatation_energy_approx_modes;
    obj.calculation_results.shear_energy_approx_modes = shear_energy_approx_modes;
    
    % check for convergence: elastic_energy and elastic_energy_approx must be similar; otherwise use a finer mesh
    
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



