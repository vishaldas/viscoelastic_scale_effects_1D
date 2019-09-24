function [ vrt,Qrt] = velrt_visco(lyr,freq,Q,omega0 )
% [ vrt,Qrt] = velrt_visco(lyr,freq,Q,omega0 )
%   Effective velocity at the RT limit
%   Input args
%
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	freq	frequency range (vector) e.g. obtained from LOGSPACE
%   Q of layered medium
%   omega0 = reference frequency for the Q of the layered medium
%   Note: include Meff in the output variables if Meff is also to be output

%   Written by Vishal Das, March 2017

% velocity dispersion equation Aki and Richards

vp=lyr(:,1); rho=lyr(:,2); thick=lyr(:,3); 
ratio = thick./sum(thick);
om = 2.*pi.*freq;

Meff = zeros(1,length(om));
vrt = zeros(1,length(om));
Qrt = zeros(1,length(om));

for i = 1:length(om)
    Qw=Q.*(1-((pi.*Q).^(-1)).*log(0.577.*om(i)./omega0));
    v=vp.*(1 +log(om(i)./omega0)./(pi.*Qw) -1i./(2.*Qw));
%     Mreal = (rho.*(real(v).^2));
%     num = sum(ratio./ (Mreal.*(1+(1./Q).^2)));
%     invqeff = qeffemt(lyr,om(i),Q,omega0);
%     Mrealeff = (num.*(1+invqeff.^2)).^(-1);
%     Meff(i) = Mrealeff.*(1+1i.*invqeff);
    denavg = sum(ratio.*rho);
%     vemt0 = sqrt(real(Meff(i))./denavg);
    vrt(i) = (sum(ratio./v)).^-1;
%     Mrt = denavg.*vrt.^2;
%     Qrt(i) = -real(Mrt)./imag(Mrt);
    wavenum = om(i)./vrt(i);
%     wavenum = real(wavenum) -1i.*(abs(imag(wavenum))); % k_imag = -alpha 
    % Qdisp = real(wavenum)./(2.*imag(wavenum)); % Carcione 2.123
    Qrt(i) = real(wavenum.^2)./(imag(wavenum.^2)); % Carcione 2.121

end

end

