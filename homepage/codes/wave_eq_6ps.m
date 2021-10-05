%% 2nd order hyperbolic PDE: the Wave Equation
% u_{tt} = c^2*u_{xx} ;  u_t(x,0)=exp(-100*(x-.7)^2)
clear all; close all; clc;

%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
r = 1.4                    % [-]   CFL (Courant number)
dt= r*dx/c;                % [s]   time sampling
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
u2 = zeros(N,1); % 'u 1 old'
uj = zeros(N,1); % 'u new'

record = [];

X = 2:N-1; % Update cells
for j=2:M
   u2(round(N/10)) =u2(round(N/10)) + fs(j)*dt^2/dx*c*2; % Inject source
   % PSEUDOSPECTRAL
   % --------> Implemented as:
   %          ifft(       -           k                             ^2   * fft(u)         )
   du = real(ifft( -([(pi/dx/N)*(0:N-1) -pi/dx+(pi/dx/N)*(0:N-1)]).^2' .* fft(u2(:),2*N)));
   du2 = real(ifft( -([(pi/dx/N)*(0:N-1) -pi/dx+(pi/dx/N)*(0:N-1)]).^2' .* fft(du(:),2*N)));
   du4 = real(ifft( -([(pi/dx/N)*(0:N-1) -pi/dx+(pi/dx/N)*(0:N-1)]).^2' .* fft(du2(:),2*N)));
   uj(X) = c^2*dt^2*du(X) + 2*c^4*dt^4/27*du2(X) + c^6*dt^6/729*du4(X)  +  2*u2(X) - u1(X);
   
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