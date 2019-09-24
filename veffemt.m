function [ veff] = veffemt(lyr,freq,Q,omega0 )
% [veff] = veffemt(  vp, rho, invq, ratio, om, omega0 )
%   Effective velocity and real part of complex modulus for low frequency limit
%   Input args
%
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	freq	frequency range (vector) e.g. obtained from LOGSPACE
%   Q of layered medium
%   omega0 = reference frequency for the Q of the layered medium
%   Note: include Meff in the output variables if Meff is also to be output

%   Written by Vishal Das, October 2016

% velocity dispersion equation Aki and Richards

vp=lyr(:,1); rho=lyr(:,2); thick=lyr(:,3); 
ratio = thick./sum(thick);
om = 2.*pi.*freq;

Meff = zeros(1,length(om));
veff = zeros(1,length(om));

for i = 1:length(om)
    Qw=Q.*(1-((pi.*Q).^(-1)).*log(0.577.*om(i)./omega0));
    v=vp.*(1 +log(om(i)./omega0)./(pi.*Qw) -1i./(2.*Qw));
    Mreal = (rho.*(real(v).^2));
    num = sum(ratio./ (Mreal.*(1+(1./Q).^2)));
    invqeff = qeffemt(lyr,om(i),Q,omega0);
    Mrealeff = (num.*(1+invqeff.^2)).^(-1);
    Meff(i) = Mrealeff.*(1+1i.*invqeff);
    denavg = sum(ratio.*rho);
    vemt0 = sqrt(real(Meff(i))./denavg);
    veff(i) = vemt0.*(2.*(1+invqeff.^2)./(sqrt(1+invqeff.^2)+1)).^(0.5);
end

end

