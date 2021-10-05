%% 2nd order hyperbolic PDE: the Wave Equation
% u_{tt} = c^2*u_{xx} ;  u_t(x,0)=exp(-100*(x-.7)^2)
clear all; close all; clc;

%% Model setup
c = 2000;                  % [m/s] wave velocity
N = 201;                   % [-]   Number of Chebyshev nodes
x = cos([0:N-1]*pi/(N-1)); % [m]   Chebyshev nodes (-1,1)
x = (x+1)*1000;            % [m]   Nodes (0,2*1000)
dx = x(1)-x(2);            % [m]   Smallest dx
dt= dx/c;                  % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
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
u2 = zeros(N,1); % 'u 1 old'
uj = zeros(N,1); % 'u new'

record = [];

X = 2:N-1; % Update cells
for j=2:M
   u2(round(N/2)) =u2(round(N/2)) + 2*fs(j)*c*dt^2/(x(round(N/2))-x(round(N/2)+1)); % Inject source

   tmp   = fastdiff(fastdiff(uj))/1000^2;
   uj(X) = c^2*dt^2*tmp(X) + 2*u2(X) - u1(X);
   
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