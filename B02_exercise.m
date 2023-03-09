%% EXERCISE #2.

%What to do:

% Design an example with 2 horizontal layers and build a plane wave 
% using a line of sources (see examples). Now change the incident angle
% of the plane wave by activating the sources in the line 
% delayed in sequence: obtain the relation between the delay, 
% incident angle, velocity. Check the Snell law. Generate a refracted wave


%% Snell's law - Critical angle

% Snell's law: sin(θ1)/V1 = sin(θ2)/V2 
% --> so we can say that sin(θi)/Vi = constant, very important for
% reflection sesismic

% We observe that the lower is V2, the smaller the transmission 
% angle and the bigger the incident angle (V2<V1).

% Critical angle: sin(θ) = V1/V2 for V1<V2. For θ = π/2 the transmitted
% wave propagated parallel to the boundary between the 2 areas. 
% There cannot exist critical angle when V1>V2.

% The model below consists of 5 different formats to illustrate different
% behaviors:
% 1: Critical angle V1=800 & V2=1200
% 2: Medium transmission angle V1=1200 & V2=800
% 3: Smaller transmission angle V1=1200 & V2=640
% 4: Bigger transmission angle V1=1200 & V2=1000
% 5: Plane wave implementation by using a line of sources

%CFL is a condition the asks if the sampling is enough for producing a wave
%equation correctly. We need to decrease dx and the better condition is
%imposing dx = dz ; or we can decrease the frequency f0 of the source.

%SHOW TWO VERSION: 1 STABLE VERSION 2 NOT STABLE VERSION, cambia il
%campionamento non più a 1 ma mettilo a 5 per avere non stabilità.
%L'indicazione che il tuotto è insatbile è che la wavefront è splittata in
%differenti wavefront e non solo 1 perchè appunto usiamo un campionamento
%troppo grande e quindi la CFL condition non è rispettata.
%Per averla stabile invece possiamo diminuire notevolmente il
%campionamento.

%%
% ----------------------------------------
% SINGLE INTERFACE WITH V2<V1
% ----------------------------------------

clear all
close all

% ----------------------------------------
% 1. Model parameters

model.x   = 0:1:1000;    % horizontal x axis sampling
model.z   = 0:1:250;     % vertical   z axis sampling

% temporary variables to compute size of velocity matrix
Nx = numel(model.x);
Nz = numel(model.z);

% example of velocity model assignement
% two layers with an interface at z_interface meters depth
z_interface=200;

for kx=1:Nx,
  x=model.x(kx);
  for kz=1:Nz,
    z=model.z(kz);
    
    if z<z_interface,
      model.vel(kz,kx)=1000;
    else
      model.vel(kz,kx)=2000;
    end
    
  end
end

% optional receivers in (recx, recz)
% the program round their position on the nearest velocity grid
model.recx  = 120:20:900;
model.recz  = model.recx*0+20;  % ... a trick to have same nr elements  of recx
model.dtrec = 0.004;
% ----------------------------------------
% 2. Source parameters

source.x    = [10:1:(model.x(end)-10)];
Nsources    = numel(source.x);

source.z    = ones(1,Nsources) * 50 ; 
source.f0   = ones(1,Nsources) * 30; 
source.t0   = ones(1,Nsources) *0.04 + (1:Nsources)*0.0008;
source.amp  = ones(1,Nsources) * 1 ;
source.type = ones(1,Nsources) * 1;    % 1: ricker, 2: sinusoidal  at f0

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
scfact = 1; % scaling factor
colour = ''; % trace colour, default is black
clip   = []; % clipping of amplitudes (if <1); default no clipping

seisplot2(recfield.data,recfield.time,[],scal,pltflg,scfact,colour,clip)
xlabel('receiver nr')


axis xy
