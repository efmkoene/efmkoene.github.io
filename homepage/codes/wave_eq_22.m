%% 2nd order hyperbolic PDE: the Wave Equation. 2nd order in time and space
% Erik Koene, ETH Z�rich, 2016
clear all; close all; clc;

%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
dt= .0015;                 % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
x = 0:dx:2000;             % [m]   cell locations
N = length(x);             % [-]   number of cells
T = 0:dt:1;                % [s]   time vector
M = length(T);             % [-]   number of time steps

% Define Ricker wavelet source
fc = 20;                   % [Hz]  central wavelet frequency
t0 = 1.5/fc;               % [s]   delay wavelet
tau=pi*(T-t0)*fc;          % [-]   Ricker wavelet argument
% fs = (1 - 2*tau.^2) .* exp( -tau.^2 );       % Original Ricker wavelet
fs = 2*pi*fc*tau .* (2*tau.^2 - 3) .* exp( -tau.^2 );    % Derivative of Ricker wavelet

% Initialize the grid
u1 = zeros(N,1); % 'u 2 old'
u2 = zeros(N,1); % 'u old'
uj = zeros(N,1); % 'u new'

record = [];

X = 2:N-1; % Update cells
for j=2:M
   u2(round(N/10)) =u2(round(N/10)) +  fs(j)*dt^2/dx*c*2; % Inject source
   % FINITE DIFFERENCE 2nd order in space, 2nd in time
   % u(t+1,x) - 2u(t,x) + u(t-1,x)       u(t,x+1) - 2u(t,x) + u(t,x-1)
   % ----------------------------- = c^2 -----------------------------
   %             dt^2                                dx^2
   %
   % -> u(t+1,x) = (c^2*dt^2/dx^2) * ( u(t,x+1) - 2u(t,x) + u(t,x-1) ) + 2u(t,x) - u(t-1,x)
   % --------> Implemented as:
   uj(X) = 2*(1-r^2)*u2(X) + r^2*(u2(X+1)+u2(X-1)) - u1(X);
   
   plot(x,uj),ylim([-1 1]),title(sprintf('time t=%0.2f',j*dt)),xlabel('offset [m]')
   hold on
   plot(x(round(3*N/5)),uj(round(3*N/5)),'*')
   hold off
   drawnow
   u1 = u2;
   u2 = uj;
   record = [record, uj(round(3*N/5))];
end

figure(2)
plot(T(1:M-1),record,'k')
title('recorded, with dispersion')
xlabel('time')