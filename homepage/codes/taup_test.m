clear all;close all;clc;

% Settings
nx = 32;
dx = 10;
xs = [0:nx-1]*dx;

nt = 80;
dt = 0.004;
ts = [0:nt-1]*dt;
% =-=-=-=-=

% Create dataset
P = zeros( nt, nx );
P(40     ,:) = P(40     ,:) + 1                ; % Flat
P(25:nt  ,:) = P(25:nt  ,:) + eye( nt-24 , nx ); % Exactly diagonal
P(10:2:nt,:) = P(10:2:nt,:) + eye( nt/2-4, nx ); % 'too' diagonal

figure(1)
imagesc(xs, ts, P )
title('Input data P(x,t)')
colorbar
% =-=-=-=-=

% Generate slowness vector
pmin = -2.5*dt/dx;
pmax =  2.5*dt/dx;
dp = 2*dt/max(xs); % From "Aliasing and tau-p" (S.I. Maroof and C.J. Gravely, Seismic 14).
ps = [pmin : dp : pmax]; % Build p vector
% =-=-=-=-=

%% Do Tau P transform backward and forward
% Tau P transform
P_tp  = taup( P   , dt, xs, ps, 1E-3,  1, 0  );
P_rec = taup( P_tp, dt, xs, ps,  0  , -1, nt );
fprintf(1,'Amplitudes P  correct up to %f %%\n', 100*max(max(abs(P-P_rec  )))/max(max(abs(P ))));
% =-=-=-=-=

% Show results
figure(2)
imagesc( ps, [0:size(P_tp,1)-1]*dt, P_tp )
title('P(\tau,p)')
colorbar

figure(3)
imagesc( xs, ts, P_rec )
title('P(x,t) reconstruction')
colorbar

figure(4)
imagesc( xs, ts, P - P_rec )
title('Error in reconstruction')
colorbar
% =-=-=-=-=