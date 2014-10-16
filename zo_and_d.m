function [zom,d] = zo_and_d (LAI,h)

% function zom_and_d calculates roughness length for momentum and zero
% plane displacement from vegetation height and LAI
%
% Date:     17 November 2008
% Author:   A. Verhoef
%           implemented into Matlab by C. van der Tol (tol@itc.nl)
%
% Source:   Verhoef, McNaughton & Jacobs (1997), HESS 1, 81-91
% 
% usage: 
%       [zom,d] = zo_and_d (LAI,h)
%
% input:
%   LAI         one sided leaf area index
%   h           vegetation height (m)
%
% output:
%   zom         roughness lenght for momentum (m)
%   d           zero plane displacement (m)
% 
% globals:
%   Cd          Averaged drag coefficient for the vegetation              
%   CR          Drag coefficient for isolated tree
%   CSSOIL      Drag coefficient for soil
%   CD1         Fitting parameter
%   Psicor      Roughness layer correction
%   kappa       Von Karman's constant

%% parameters
global CR CSSOIL CD1 Psicor kappa

%% calculations
sq      = sqrt(CD1*LAI/2);
G1      = max(3.3, (CSSOIL + CR*LAI/2).^(-0.5));
d       = (LAI>1E-7 & h>1E-7).*h.*(1-(1-exp(-sq))./sq);          % Eq 12 in Verhoef et al (1997)
zom     = (h-d).*exp(-kappa*G1 + Psicor);


