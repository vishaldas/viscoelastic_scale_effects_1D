function [ invqeff,veff ] = viscoemt( lyr,Q)
%[ invqeff,veff ] = viscoemt( lyr,Q)
%   Effective invq and velocity for low frequency limit based on just Q and layer
%   properties for a particular frequency. 
%   Input args
%
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%   Q of layered medium

%   Written by Vishal Das, November 2017



vp=lyr(:,1); rho=lyr(:,2); thick=lyr(:,3);
ratio = thick./sum(thick);

Mreal = real(rho.*vp.^2);
num = sum((ratio.*(1./Q))./(Mreal.*(1+(1./Q).^2)));
den = sum(ratio./(Mreal.*(1+(1./Q).^2)));
invqeff = num/den;

num1 = sum(ratio./ (Mreal.*(1+(1./Q).^2)));
Mrealeff = (num1.*(1+invqeff.^2)).^(-1);
Meff = Mrealeff.*(1+1i.*invqeff);
denavg = sum(ratio.*rho);
vemt0 = sqrt(real(Meff)./denavg);
veff = vemt0.*(2.*(1+invqeff.^2)./(sqrt(1+invqeff.^2)+1)).^(0.5);


end

