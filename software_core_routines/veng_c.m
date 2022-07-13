
function [ff,phi_ff] = veng_c(a,b,n,params,T)

    % Calculates thermoelastic dissipation for a two-dimensional structure, 
    % according to Vengallatore's model,
    % in a sample WITH COATING.
    % a, b = geometrical parameters as defined in Vengallatore's model
    % n = number of terms in the series to be calculated
    % T = temperature [K]
    % params = structure containing the material properties (see comments in this file for details).
    % This implementation is inspired by the work of M. Lorenzini.
    %
    % a, b, T - double.
    % n - integer.
    % params - struct.
    % ff, phi_ff - array 1x5000 double. 

    Cs = params.Cs;             % heat_capacity = specific_heat_capacity*density (of Substrate)
    ks = params.ks;             % thermal_conductivity (of Substrate)
    alphas = params.alphas;     % linear_thermal_expansion (of Substrate)
    Es = params.Es;             % Young's modulus (of Substrate)

    Cc = params.Cc;             % heat_capacity = specific_heat_capacity*density (of Coating)
    kc = params.kc;             % thermal_conductivity (of Coating) 
    alphac = params.alphac;     % linear_thermal_expansion (of Coating) 
    Ec = params.Ec;             % Young's modulus (of Coating)  
    
    %%%%

    taus = Cs/ks*a^2;
    tauc = Cc/kc*(b-a)^2;

    %%%%

    % Find equation roots:
    % first searching for minima of fdiff (with a sampling step independent of n) 
    % then checking and refining with vpasolve (using the minima as initial guesses).
    % There should be about n roots in the interval [0,n*pi/sqrt(taus)]

    K = sqrt(ks*Cs/kc/Cc);
    fdiff = @(xx) tan(xx*sqrt(taus)) - K/tan(xx*sqrt(tauc));

    syms xsym
    eqn =  fdiff(xsym) == 0;

    step = (pi/sqrt(taus))/10000000;
    beta_n = [];
    posxx = 0;
    while length(beta_n) < n
        if abs(fdiff(posxx)) > abs(fdiff(posxx+step)) && abs(fdiff(posxx+step)) < abs(fdiff(posxx+2*step))
            beta_n(end+1) = vpasolve(eqn,xsym,posxx+step);
        end
        posxx = posxx + step;
    end

    %%%%
    
    gamma_n = beta_n*sqrt(taus);
    eta_n = beta_n*b/(b-a)*sqrt(tauc);

    A_n = - K*cos(gamma_n).*cos(eta_n)./sin(eta_n*(a/b-1));
    B_n = - K*cos(gamma_n).*sin(eta_n)./sin(eta_n*(a/b-1));

    I1_n = 1./gamma_n.^2.*sin(gamma_n) - 1./gamma_n.*cos(gamma_n);
    I2_n = 1./eta_n.*(...
        A_n .* (1./eta_n.*(cos(eta_n)-cos(eta_n*a/b)) + sin(eta_n) - a/b*sin(eta_n*a/b)) + ...
        B_n .* (1./eta_n.*(sin(eta_n)-sin(eta_n*a/b)) - cos(eta_n) + a/b*cos(eta_n*a/b)) ...
        );
    I3_n = 1/2 - 1./(4*gamma_n).*sin(2*gamma_n);
    I4_n = 1/2*(1-a/b) * (A_n.^2 + B_n.^2) + ...
        1./(4*eta_n) .* (A_n.^2 - B_n.^2) .* (sin(2*eta_n)-sin(2*eta_n*a/b)) - ...
        1./(2*eta_n) .* A_n.*B_n .* (cos(2*eta_n)-cos(2*eta_n*a/b));

    Q_n = ( (a^2/b^2*I1_n + Ec*alphac/(Es*alphas)*I2_n).^2 ) ./ ...
          ( a/b*I3_n + Cc/Cs*I4_n );

    %%%%

    psi0s = alphas^2*T*Es/Cs;               % Zener's modulus (of Substrate)
    omegapeak_n = ks/(Cs*a^2)*gamma_n.^2;    % 2*pi*f value of max position 
    
    %%%%%
    
    ff = logspace(1,6,5000); 
    oo = 2*pi*ff;

    phi_ff = zeros(1,length(ff));
    for mm = 1:n
        phi_ff = phi_ff + ...
            psi0s / (1/3*(a^3/b^3 + Ec/Es*(1-a^3/b^3))) * ...
            ( (oo*omegapeak_n(mm)) ./ (oo.^2 + omegapeak_n(mm)^2) * Q_n(mm) );
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




