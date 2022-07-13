
function [v1, v2] = find_orthogonal_vectors(v3)

    % Given in input the coordinates of one vecor
    % (not necessarily normalised in legth), expressed in an orthonormal basis,
    % returns the coordinates of two normalised vectors forming a right-handed triad
    % with the one given as input (which corresponds to the third vector).
    % The two perpendicular directions are aribtrarily chosen.
    %
    % v3 - array 1x3 or 3x1 double.
    % v1, v2 - arrays 1x3 double.

    v3 = reshape(v3,[1,3]);

    v3 = v3/norm(v3);

    null_space_basis = null(v3(:)');

    v1 = null_space_basis(:,1);
    v2 = null_space_basis(:,2);

    v1 = reshape(v1,[1,3]);
    v2 = reshape(v2,[1,3]);

    % check right-handedness with: 
    % disp([num2str(cross(v1,v2)) ' must be = ' num2str(v3)]);
    % disp([num2str(cross(v2,v3)) ' must be = ' num2str(v1)]);
    % disp([num2str(cross(v3,v1)) ' must be = ' num2str(v2)]);
    
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





