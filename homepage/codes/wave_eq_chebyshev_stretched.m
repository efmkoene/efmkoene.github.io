%% 2nd order hyperbolic PDE: the Wave Equation
% u_{tt} = c^2*u_{xx} ;  u_t(x,0)=exp(-100*(x-.7)^2)
clear all; close all; clc;

%% Model setup
c = 2000;                  % [m/s] wave velocity
N = 301;                   % [-]   Number of Chebyshev nodes
e = cos([0:N-1]*pi/(N-1)); % [m]   Chebyshev nodes (-1,1)
a = sech(log(abs(log(eps)))/N);  % [-]   Stretch-regularization choice [0,1]
x = 1000*(1 - asin(a*e)/asin(a) );  % [m] Stretched grid (regularized)
dude = asin(a)/a/1000 * sqrt(1-(a*e).^2);
dx = x(2)-x(1);            % [m]   Smallest dx
dt= 1*dx/c;                  % [s]   time sampling
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

X = 1:N; % Update cells
for j=2:M
   u2(round(N/2)-20) =u2(round(N/2)-20) + 1.95*fs(j)*c*dt^2/(x(round(N/2)-19)-x(round(N/2)-20)); % Inject source

   tmp   = fastdiff(uj); tmp=tmp.*dude'; tmp( [1 N] ) = 0; % boundary condition: p'=0 at the boundaries!
   tmp   = fastdiff(tmp).*dude';
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