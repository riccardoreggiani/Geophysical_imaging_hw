% -------------------------------------------------------------------------
%       Acoustic wave equation finite difference simulator
% -------------------------------------------------------------------------

% --------------------------------------------------
% In a homogeneous medium, compute the energy of traceS at incresaing 
% distance from the source
% spherical divergence
% Plot the result energy vs distance

clear all


% ----------------------------------------
% 1. Model parameters

model.x   = 0:1:1000;    % horizontal x axis sampling
model.z   = 0:1:450;     % vertical   z axis sampling

% temporary variables to compute size of velocity matrix
Nx = numel(model.x);
Nz = numel(model.z);

% example of velocity model assignement
% two layers with an interface at z_interface meters depth
z_interface=200;
model.vel=zeros(Nz,Nx);      % initialize matrix


% Homogeneous medium


for kx=1:Nx,
  x=model.x(kx);
  for kz=1:Nz,
    z=model.z(kz);
    model.vel(kz,kx)=1300;  
  end 
    end


% optional receivers in (recx, recz)
% the program round their position on the nearest velocity grid
model.recx  = 120:20:900;
model.recz  = model.recx*0+20;  % ... a trick to have same nr elements  of recx
model.dtrec = 0.004;
% ----------------------------------------
% 2. Source parameters

source.x    = 100;
Nsources    = numel(source.x);

source.z    = 20 ; 
source.f0   = ones(1,Nsources) * 30; 
source.t0   = ones(1,Nsources) * 0.04;
for n=2:Nsources
    source.t0(n)   = source.t0(n) + 0.0001*n;
end
source.amp  = ones(1,Nsources) * 1 ;
source.type = ones(1,Nsources) * 1;    % 1: ricker, 2: sinusoidal  at f0


% ----------------------------------------
% 3. Simulation and graphic parameters in structure simul

simul.borderAlg=1;
simul.timeMax=1;

simul.printRatio=10;
simul.higVal=.05;
simul.lowVal=0.02;
simul.bkgVel=1;

simul.cmap='gray';   % gray, cool, hot, parula, hsv

% ----------------------------------------
% 4. Program call

recfield=acu2Dpro(model,source,simul);



% Plot receivers traces

figure
scal   = 2;  % 1 for global max, 0 for global ave, 2 for trace max
pltflg = 0;  % 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
             % 2 plot wiggle traces only, 3 imagesc gray, 4 pcolor gray
scfact = 5; % scaling factor
colour = ''; % trace colour, default is black
clip   = []; % clipping of amplitudes (if <1); default no clipping

seisplot2(recfield.data,recfield.time,[],scal,pltflg,scfact,colour,clip)
xlabel('receiver nr')
axis xy


%energy traces
energy = recfield.data .^2;
figure
scal   = 2;  % 1 for global max, 0 for global ave, 2 for trace max
pltflg = 0;  % 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
             % 2 plot wiggle traces only, 3 imagesc gray, 4 pcolor gray
scfact = 5; % scaling factor
colour = ''; % trace colour, default is black
clip   = []; % clipping of amplitudes (if <1); default no clipping

seisplot2(energy,recfield.time,model.recx-source.x,scal,pltflg,scfact,colour,clip)
xlabel('distance [m]')
axis xy

%enrgy attenuation
dist = model.recx-source.x;
massimi = max(energy, [], 1);
figure
plot(dist,massimi)
xlabel('Distance [m]')


