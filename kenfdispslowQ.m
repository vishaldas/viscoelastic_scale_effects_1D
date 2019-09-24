function [ freq, veldisp, vel_rt, vel_emt,Qdisp ] = kenfdispslowQ( lyr,freq,Q,omega0 )
%function [ freq, veldisp, vel_rt, vel_emt,Qdisp ] = kenfdispslowQ( lyr,freq,Q,omega0 )
%[ freq, veldisp, vel_rt, vel_emt ] = kenfdispslowQ( lyr,freq,Q,omega0 )
%   Kennett-Frazer dispersion curve (scattering + intrinsic) using slowness 
%   Velocity dispersion for 1-D layered media, normal incidence plane wave
%   using the Kennett-Frazer algorithm. This is the scattering and
%   intrinsic, viscoelastic dispersion.
%
%   Inputs
%	LYR	[velocity, density, thickness] of layered medium
%		nx3 matrix, n=number of layers.
%	freq	frequency range (vector) e.g. obtained from LOGSPACE
%   Q of layered medium
%   omega0 = reference frequency for the Q of the layered medium
%   Outputs
%	VELDISP	dispersion velocity (vector of length = length of F).
%   VEL_RT  velocity at ray theory limit of viscoelastic medium
%   VEL_EMT velocity at EMT limit of viscoelastic medium
%   QDISP calculated Q from the dispersion velocity 
%
%   Written by Vishal Das, January 2017


vel = lyr(:,1); rho = lyr(:,2); thick = lyr(:,3);
n = length(vel);

om = 2.*pi.*freq;
[~,om2] = size(om);

% Calculating the RT limit
time_rt = sum(repmat(thick,1,om2)./repmat(vel,1,om2));
vel_rt = sum(repmat(thick,1,om2))./time_rt;

% Calculating the EMT limit
vel_emt = veffemt(lyr,freq(1),Q,omega0);

f = thick./sum(thick);
f = repmat(f,1,length(freq));

[v1,~] = size(vel);
om = repmat(om,v1,1);
vel  = repmat(vel,1,om2);
rho = repmat(rho,1,om2); 
Q=repmat(Q,1,om2); 

Qw=Q.*(ones(v1,om2)-((pi.*Q).^(-1)).*log(0.577.*om./omega0));
Qw(isinf(Qw)) = mean(Qw(~isinf(Qw)));
% Modified velocity, from Aki and Richards equation 5.94
vel=vel.*(ones(v1,om2) +log(om./omega0)./(pi*Qw) -1i./(2*Qw));
vel(isinf(vel)) = mean(vel(~isinf(vel)));

deno = ones(size(vel)); rd = ones(size(vel)); td = ones(size(vel));
deno(2:n,:) = rho(2:n,:).*(vel(2:n,:)) + rho(1:n-1,:).*(vel(1:n-1,:));
rd(2:n,:)   = (rho(1:n-1,:).*(vel(1:n-1,:)) - rho(2:n,:).*(vel(2:n,:)))./deno(2:n,:);
td(2:n,:)   = 2*sqrt(rho(2:n,:).*(vel(2:n,:)).*rho(1:n-1,:).*(vel(1:n-1,:)))./deno(2:n,:);

% Considering free surface at the top
rd(1,:) = -1*rd(1,:);
td(1,:) =  1*td(1,:);

ru=-rd; tu=td;

% Slowness calculation from velocity
S = 1./vel;
% Considering positive part of the wave propagaion
% S = real(S)+1i.*abs(imag(S));

Ru =0; Rd=0; Td=1; Tu=1;
xsum= 0; 

for k = n:-1:1
    theta = exp(1i.*(thick(k).*S(k,:)).*om(k,:));
    xsum = xsum + log(td(k,:).*((1-Rd.*theta.^2.*ru(k,:)).^-1));
    % xsum = xsum + log(tu(k).*((1-Rd.*theta.^2.*ru(k)).^-1));
    % Ru = (-r + Ru.*theta.^2)./(1-r*Ru.*theta.^2); % kenfdisp code
    Ru = -rd(k,:)+td(k,:).*(1-Ru.*theta.^2.*rd(k,:)).^(-1).*Ru.*theta.^2.*tu(k,:); % Frazer 1994
    % Ru = ru(k)+td(k).*(1-Ru.*theta.^2.*rd(k)).^(-1).*Ru.*theta.^2.*tu(k); % Frazer 1994
    reverb = 1./(1-Rd.*theta.^2.*ru(k,:));
    Rd=rd(k,:)+tu(k,:).*theta.*Rd.*theta.*reverb.*td(k,:);
    Td=Td.*theta.*reverb.*td(k,:);
end

Sin = sum(S(1:end,:).*f(1:end,:));
Sst = (1./(1i.*om(1,:).*sum(thick(1:end)))).*xsum;
Stotal = (Sst) + Sin;
veldisp = 1./real(Stotal);

wavenum = Stotal.*2.*pi.*freq;
% wavenum = real(wavenum) -1i.*(abs(imag(wavenum))); % k_imag = -alpha 
% Qdisp = real(wavenum)./(2.*imag(wavenum)); % Carcione 2.123
Qdisp = real(wavenum.^2)./(imag(wavenum.^2)); % Carcione 2.121

