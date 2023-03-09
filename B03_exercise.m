% -------------------------------------------------------------------------
%       Acoustic wave equation finite diference simulator
% -------------------------------------------------------------------------

% ----------------------------------------
% Design an example with some horizontal layers, and collect a common 
% shot gather with a line of receivers in the first layer.
% Verify that the reflecton traveltime curve for each interface is a
% hyperbola of equation
% t_k(x)=sqrt(t_0k^2 + (x/V_rms_k))^2)
% where 
% k : interface nr
% t_0k  : two way traveltime for vertical propagation from recording level
%         to interface k
% V_rms_k=sqrt( sum_for_n=1_to_k(V_n^2*TWT_n)/
%               sum_for_n=1_to_k(TWT_n)
%              )
% V_n : velocity of layer n
% TWT_n : two way traveltime for vertical propagation *within* layer n
%
% trace the traveltime curves over the shot gather


clear all

% 1. Model parameters

model.x   = 0:1:800;     % horizontal x axis sampling
model.z   = 0:1:500;     % vertical   z axis sampling

% temporary variables to compute size of velocity matrix
Nx = numel(model.x);
Nz = numel(model.z);

% velocity model assignement and layers

for kx=1:Nx
    x=model.x(kx);
    for kz=1:Nz
        z=model.z(kz);
        if z<200
            model.vel(kz,kx)=1000;
        elseif 200<=z && z<300
            model.vel(kz,kx)=1500;
        elseif 300<=z && z<400
            model.vel(kz,kx)=2000;
        elseif z>=400
            model.vel(kz,kx)=2500;
        end
    end
end

% ----------------------------------------
% 2. Source parameters


source.x    = [100];
source.z    = [100]; 
source.f0   = [25];
source.t0   = [0.04];
source.amp  = [5];
source.type = [1];    % 1: ricker, 2: sinusoidal  at f0

% optional receivers in (recx, recz)
% the program round their position on the nearest velocity grid

model.recx  = [200:10:(model.x(end)-100)];
Nreceivers  = numel(model.recx);

model.recz  = ones(1,Nreceivers) * 100;
model.dtrec = 0.004;

% ----------------------------------------
% 3. Simulation and graphic parameters in structure simul

simul.borderAlg=1;
simul.timeMax=1;

simul.printRatio=10;
simul.higVal=.6;
simul.lowVal=.1;
simul.bkgVel=1;

simul.cmap='gray';   % gray, cool, hot, parula, hsv


% ----------------------------------------
% 4. Calculations

V_n = [1000,1500,2000,2500];
z=[100,200,300,400];
depth=100;
TWT_n = 2*depth./V_n;
k=length(V_n)-1;
V_rms_k = zeros(1,k);
for i=1:k
    sumation1=0;
    sumation2=0;
    for j=1:i
        sumation1 = sumation1 + (V_n(j))^2*TWT_n(j);
        sumation2 = sumation2 + TWT_n(j);
    end 
    V_rms_k(i) = sqrt(sumation1/sumation2);
end

t0_k = size(V_rms_k);
for i=1:k
    t0_k(i) = 2*z(i)/V_rms_k(i);
end

t_k = zeros(Nx,k);
for i=1:k
    for j=1:Nx
        t_k(j,i)=source.t0 + sqrt((t0_k(i))^2+(j/V_rms_k(i))^2);
    end
end



% ----------------------------------------
% 5. Program call

recfield=acu2Dpro(model,source,simul);

% Plot receivers traces

recfield.offset = (model.recx - source.x);

figure
scal   = 1;  % 1 for global max, 0 for global ave, 2 for trace max
pltflg = 0;  % 1 plot only filled peaks, 0 plot wiggle traces and filled peaks,
             % 2 plot wiggle traces only, 3 imagesc gray, 4 pcolor gray
scfact = 6;  % scaling factor
colour = ''; % trace colour, default is black
clip   = []; % clipping of amplitudes (if <1); default no clipping

seisplot2(recfield.data,recfield.time,recfield.offset,scal,pltflg,scfact,colour,clip)
xlabel('Offset[m]')
hold on
plot(t_k,'r')