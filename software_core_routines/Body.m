
classdef Body 
    
    % Class containing information on both the material and geometric properties; 
    % instantiated with the object name "substrate" or "coating".
    
    properties
       
        thickness
        diameter
        material
        
    end
    
    methods
        
        function obj = Body()
            
            % constructor function, automatically called at instantiation
            
            obj.material = Material();
            
        end
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

