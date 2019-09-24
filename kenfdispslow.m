function [ freq, veldisp, vel_rt, vel_emt,Qdisp ] = kenfdispslow( lyr,freq )
%function [ freq, veldisp, vel_rt, vel_emt ] = kenfdispslow( lyr,freq )
%[ freq, veldisp, vel_rt, vel_emt ] = kenfdispslow( lyr,freq )
%   Kennett-Frazer dispersion curve using slowness 
%   Velocity dispersion for 1-D layered media, normal incidence plane wave
%   using the Kennett-Frazer algorithm. This is the scattering dispersion,
%   not intrinsic, viscoelastic dispersion.
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	freq	frequency range (vector) e.g. obtained from LOGSPACE
%   Outputs
%	VELDISP	dispersion velocity (vector of length = length of F).
%   VEL_RT  velocity at ray theory limit
%   VEL_EMT velocity at EMT limit
%
%   Written by Vishal Das, January 2017


vel = lyr(:,1); rho = lyr(:,2); thick = lyr(:,3); 
w = 2.*pi.*freq;
f = thick./sum(thick);
n = length(vel);

% Calculating the RT limit
time_rt = sum(thick./vel);
vel_rt = sum(thick)./time_rt;

% Calculating the EMT limit
m = rho.*vel.^2;
den_avg = sum(f.*rho);
memt = (sum(f./m)).^-1;
vel_emt = sqrt(memt./den_avg);
time_emt = sum(thick)./vel_emt;

% Calculating the reflection and transmission coeffs
j=2:n;
deno=rho(j).*vel(j)+rho(j-1).*vel(j-1);
rd=(rho(j-1).*vel(j-1)-rho(j).*vel(j))./deno;
td=2*sqrt(rho(j).*vel(j).*rho(j-1).*vel(j-1))./deno;
rd=[-1;rd]; td=[1;td];
% rd(end) = 0; td(end) = 1;
ru=-rd; tu=td;

S = 1./vel;

Ru =0; Rd=0; Td=1; Tu=1;
xsum= 0; 

for k = n:-1:1
    theta = exp(1i.*(thick(k).*S(k)).*w);
    xsum = xsum + log(td(k).*((1-Rd.*theta.^2.*ru(k)).^-1));
    % xsum = xsum + log(tu(k).*((1-Rd.*theta.^2.*ru(k)).^-1));
    % Ru = (-r + Ru.*theta.^2)./(1-r*Ru.*theta.^2); % kenfdisp code
    Ru = -rd(k)+td(k).*(1-Ru.*theta.^2.*rd(k)).^(-1).*Ru.*theta.^2.*tu(k); % Frazer 1994
    % Ru = ru(k)+td(k).*(1-Ru.*theta.^2.*rd(k)).^(-1).*Ru.*theta.^2.*tu(k); % Frazer 1994
    reverb = 1./(1-Rd.*theta.^2.*ru(k));
    Rd=rd(k)+tu(k).*theta.*Rd.*theta.*reverb.*td(k);
    Td=Td.*theta.*reverb.*td(k);
end

Sst = (1./(1i.*w.*sum(thick(1:end)))).*xsum;
Stotal = (Sst) + (1./vel_rt);

veldisp = 1./real(Stotal);

wavenum = Stotal.*2.*pi.*freq;
% wavenum = real(wavenum) -1i.*(abs(imag(wavenum))); % k_imag = -alpha 
% Qdisp = real(wavenum)./(2.*imag(wavenum)); % Carcione 2.123
Qdisp = real(wavenum.^2)./(imag(wavenum.^2)); % Carcione 2.121
