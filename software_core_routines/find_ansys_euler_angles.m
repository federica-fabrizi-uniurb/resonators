
function out = find_ansys_euler_angles(a,b,c)

    % Given the coordinates of three orthogonal vectors a, b, c forming a right-handed triad
    % (not necessarily normalised in legth), expressed in an orthonormal basis,
    % returns the Euler angles (in radians)
    % bringing the triad into the basis vectors x, y, z:
    % a into x = [1 0 0]; b into y = [0 1 0]; c into z = [0 0 1].
    % This is used to specify the orientation of a single crystal material with
    % respect to the geometrical shape of the resonator.
    % The Euler rotations correspond to an intrinsic rotation about z, then
    % about x' (x axis rotated by the first rotation), then 
    % about y'' (y axis rotated by the first and the second rotations),
    % as requested by Ansys for the MAPDL input file "input_file_edited.dat".
    % See cooments in this file for details and file Material.m for conventions on rotations.
    %
    % a, b, c - arrays 1x3 or 3x1 double.
    % out - array 1x3 double.

    % What we want:
    % After the rotations (consider them Active rotations, where the basis vectors are defined by coordinates and 
    % are objects to be rotated),
    % we want:
    % R*[c(1) c(2) c(3)]' = R*c' = [0 0 1]';
    % and also
    % R*[a(1) a(2) a(3)]' = R*a' = [1 0 0]';
    
    % R*[b(1) b(2) b(3)]' = R*b' = [0 1 0]';

    % Therefore:

    a = reshape(a,[3,1]);
    b = reshape(b,[3,1]);
    c = reshape(c,[3,1]);

    a = a/norm(a);
    b = b/norm(b);
    c = c/norm(c);

    R = inv([a(1) b(1) c(1); a(2) b(2) c(2); a(3) b(3) c(3)]);

    % Now we want to express R as a sequence of 3 Euler rotations.

    % Ansys defines the new reference system starting from the global system
    % and performing 3 INTRINSIC rotations
    % 1st) about z   with phi_z > 0 for right-handed rotation of the new system
    % 2nd) about x'  with phi_x > 0 for right-handed rotation of the new system
    % 3rd) about y'' with phi_y > 0 for right-handed rotation of the new system

    % R = M_z(phi_z) * M_x(phi_x) * M_y(phi_y)
    % A sequence of 3 intrinsic rorations, in order:
    %    1st abiut z   by phi_z
    %    2nd about x'  by phi_x
    %    3rd about y'' by phi_y

    % However, MatLab does not support this order:
    % "The given Euler sequence ZXY is not supported. Try one of the supported sequences: ZYX, ZYZ, XYZ."

    % Then, we RE-NAME THE AXES so that 
    % X_ANSYS = y_matlab
    % Y_ANSYS = z_matlab
    % Z_ANSYS = x_matlab

    % Which means we change R to Rprime:

    Rprime = [[R(3,3) R(3,1) R(3,2)];[R(1,3) R(1,1) R(1,2)];[R(2,3) R(2,1) R(2,2)]];

    eul_matlab = rotm2eul(Rprime,'XYZ');

    % as a check:

    psi_x = eul_matlab(1);
    psi_y = eul_matlab(2);
    psi_z = eul_matlab(3);

    Mp_x = [[1   0             0          ];
            [0   +cos(psi_x)   -sin(psi_x)];
            [0   +sin(psi_x)   +cos(psi_x)]];

    Mp_y = [[+cos(psi_y)   0   +sin(psi_y)];
            [0             1   0          ];
            [-sin(psi_y)   0   +cos(psi_y)]];   

    Mp_z = [[+cos(psi_z)   -sin(psi_z)   0];       
            [+sin(psi_z)   +cos(psi_z)   0];
            [0             0             1]];    

    Rprime_check =  Mp_x * Mp_y * Mp_z;   

    % must be Rprime_check = Rprime.

    % We can interpret this sequence as the 3 rotations we wanted, in the
    % reference frame with axes not re-named:

    % BACK IN ANSYS AXES:

    % rotation about x_matlab = rotation about Z_ANSYS
    % rotation about y_matlab = rotation about X_ANSYS
    % rotation about z_matlab = rotation about Y_ANSYS

    % so the three angles are what we need.

    % as a check:

    phi_z = eul_matlab(1);
    phi_x = eul_matlab(2);
    phi_y = eul_matlab(3);

    M_x = [[1   0             0          ];
           [0   +cos(phi_x)   -sin(phi_x)];
           [0   +sin(phi_x)   +cos(phi_x)]];

    M_y = [[+cos(phi_y)   0   +sin(phi_y)];
           [0             1   0          ];
           [-sin(phi_y)   0   +cos(phi_y)]];   

    M_z = [[+cos(phi_z)   -sin(phi_z)   0];       
           [+sin(phi_z)   +cos(phi_z)   0];
           [0             0             1]];    

    R_check =  M_z * M_x * M_y;

    % must be R_check = R.

    out = eul_matlab;     
    
end

% See also:

% 1) "Here we present the results for the two most commonly used
%     conventions: ZXZ for proper Euler angles and ZYX for Tait–Bryan. 
%     Notice that any other convention can be obtained just changing the name of the axes."
%    https://en.wikipedia.org/wiki/Euler_angles#Matrix_orientation   

% 2) https://web.archive.org/web/20220628132852/https://stackoverflow.com/questions/6770845/how-to-convert-euler-xyz-to-euler-zxy-using-javascript
%     archived from
%    https://stackoverflow.com/questions/6770845/how-to-convert-euler-xyz-to-euler-zxy-using-javascript




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



