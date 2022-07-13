
function [ff,phi_ff] = veng_s(a,n,params,T)

    % Calculates thermoelastic dissipation for a two-dimensional structure, 
    % according to Vengallatore's model,
    % in a sample made of BARE SUBSTRATE (NO COATING).
    % a = geometrical parameter as defined in Vengallatore's model (b assumed to be = 0) 
    % n = number of terms in the series to be calculated
    % T = temperature [K]    
    % params = structure containing the material properties (see comments in this file for details).
    % This implementation is inspired by the work of M. Lorenzini.
    %
    % a, T - double.
    % n - integer.
    % params - struct.
    % ff, phi_ff - array 1x5000 double.     

    Cs = params.Cs;             % heat_capacity = specific_heat_capacity*density (of Substrate)
    ks = params.ks;             % thermal_conductivity (of Substrate)
    alphas = params.alphas;     % linear_thermal_expansion (of Substrate)
    Es = params.Es;             % Young's modulus (of Substrate)
    
    %%%%

    gamma_n = [];
    for mm = 1:n
        gamma_n(end+1) = (2*mm-1)*pi/2;
    end

    %%%%

    psi0s = alphas^2*T*Es/Cs;               % Zener's modulus (of Substrate)
    omegapeak_n = ks/(Cs*a^2)*gamma_n.^2;    % 2*pi*f value of max position 

    %%%%%
    
    ff = logspace(1,6,5000); 
    oo = 2*pi*ff;

    phi_ff = zeros(1,length(ff));
    for mm = 1:n
        phi_ff = phi_ff + ...
            psi0s * 6 * ...
            oo * omegapeak_n(mm) ./ (oo.^2 + omegapeak_n(mm)^2) * ...
            1./gamma_n(mm)^4 ;
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




