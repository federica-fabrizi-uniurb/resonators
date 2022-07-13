
function volumechange = get_ansys_volumechange(obj,mode_number)

    % Loads from the .txt files stored in output_files, for a given mode, the
    % element-by-element results for the (undistorted) volumes and the
    % strains, as calculated by Ansys,
    % then uses them to calculate the difference between distorted and
    % undistorted volumes (see comments in this file for details).
    %
    % obj - instance of class Disc, DoublyCoatedDisc, CurvedDisc or CantileverFiber (or any other subclass of Resonator).
    % mode_number - integer.
    % volumechange - array Nx1 double (N = number of elements).    
    
    % trace of strain tensor = delta V / V_undist = (V_dist - V_undist) / V_undist
    % The quantity returned here is: 
    % volumechange = V_dist - V_undist
    % from which one can find:
    % V_dist = volumechange + V_undist = volumechange + volume
    % (valid for small distortions)
        
    [strainXX, strainYY, strainZZ, strainXY, strainYZ, strainXZ] = get_ansys_strain(obj,mode_number);
    volume = get_ansys_volume(obj,mode_number);
    
    volumechange = (strainXX + strainYY + strainZZ).*volume;
    
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



