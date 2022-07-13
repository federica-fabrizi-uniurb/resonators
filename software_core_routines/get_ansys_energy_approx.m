
function [elastic_energy_approx,dilatation_energy_approx,shear_energy_approx] = get_ansys_energy_approx(obj,mode_number)

    % Loads from the .txt files stored in output_files, for a given mode, the
    % element-by-element results for stress, strain and volume, and uses
    % them to compute the element-by-element approximated values of the elastic energy,
    % separated into dilatation and shear contributions.
    % For each element, it must hold that:
    % elastic_energy_approx = dilatation_energy_approx + shear_energy_approx.
    % It is approximated in that it does not perform the full element integration
    % but uses the values at the centroid of the element.
    % The values returned for elastic_energy_approx can be comapred with those returned by function
    % get_ansys_sene, to judge the quality of the mesh, and consequently
    % the reliability of the values returned for dilatation_energy_approx and shear_energy_approx.
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    % mode_number - integer.    
    % elastic_energy_approx, dilatation_energy_approx, shear_energy_approx - arrays Nx1 double (N = number of elements).

    [strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = obj.get_element_strain(mode_number);
    [stressXX, stressYY, stressZZ, stressXY, stressYZ, stressXZ] = obj.get_element_stress(mode_number);
    volume = obj.get_element_volume(mode_number);
    
    elastic_energy_approx = 1/2 * (strainXX.*stressXX + strainYY.*stressYY + strainZZ.*stressZZ ...
                               + strainXY.*stressXY + strainYZ.*stressYZ + strainXZ.*stressXZ) .* volume;
                                             
    dilatation_energy_approx = 1/6 * (strainXX.*stressXX + strainYY.*stressYY + strainZZ.*stressZZ ...
                               + strainXX.*stressYY + strainXX.*stressZZ ...
                               + strainYY.*stressXX + strainYY.*stressZZ ...
                               + strainZZ.*stressXX + strainZZ.*stressYY) .* volume;
                                                 
    shear_energy_approx = 1/2 * ( 2/3*(strainXX.*stressXX + strainYY.*stressYY + strainZZ.*stressZZ) + ...
                                   strainXY.*stressXY + strainYZ.*stressYZ + strainXZ.*stressXZ + ...
                                   - 1/3*( strainXX.*stressYY + strainXX.*stressZZ ...
                                         + strainYY.*stressXX + strainYY.*stressZZ ...
                                         + strainZZ.*stressXX + strainZZ.*stressYY) ) .* volume;                     
                           
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


