function [ invqeff ] = qeffemt( lyr,freq,Q,omega0)
%[ invqeff ] = qeffemt( lyr,freq,Q,omega0)
%   Effective invq for low frequency limit for standard linear solid model
%   Input args
%
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	freq	frequency range (vector) e.g. obtained from LOGSPACE
%   Q of layered medium
%   omega0 = reference frequency for the Q of the layered medium

%   Written by Vishal Das, September 2016

% velocity dispersion equation Aki and Richards

vp=lyr(:,1); rho=lyr(:,2); thick=lyr(:,3);
ratio = thick./sum(thick);
om = 2.*pi.*freq;
invqeff = zeros(1,length(om));

for i = 1: length(om)
    Qw=Q.*(1-((pi.*Q).^(-1)).*log(0.577.*om(i)./omega0));
    v=vp.*(1 +log(om(i)./omega0)./(pi.*Qw) -1i./(2.*Qw));
    Mreal = real(rho.*v.^2);
    num = sum((ratio.*(1./Q))./(Mreal.*(1+(1./Q).^2)));
    den = sum(ratio./(Mreal.*(1+(1./Q).^2)));
    invqeff(i) = num/den;
end
end

