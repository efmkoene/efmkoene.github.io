%% 1st order hyperbolic PDE: the Wave Equation
clear all; close all; clc;

% Process gif with Bash command: gifsicle -b -O3 -d3 --lossy=80 --colors 40 --crop 25,0+-25x-0 staggered_grid.gif


%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
dt= .0015;                 % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
x = 0:dx:2000;             % [m]   cell locations
N = length(x);             % [-]   number of cells
T = 0:dt:1;                % [s]   time vector
M = length(T);             % [-]   number of time steps
rho=1500;

% Define Ricker wavelet source
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
f      = -ricker(30 , T-0.06            ) / 2; % Monopole source
g      =  ricker(30 , T-0.06+(dx/c-dt)/2) / 2; % Dipole source
f=FTDT(f);
g=FTDT(g);


% Initialize the pressure and particle velocity grid
[p,v] = deal( zeros(N,1) );

% Seismogram init
record = zeros(M,1);

% GIF INIT
h=figure(1); set(gcf,'Position',[2700 200 562 225]); axis tight;
filename='../images/staggered_grid.gif';
gif = 1;

% Optimal FD coefficients, made with: co=FD_coeffs(2,1,5e-4)
co = FD_coeffs(2,1,1e-16);

XX = 3:N-2; % Update cells
for j=2:M
    % MONOPOLE SOURCE 
    p(round(N/10)) =p(round(N/10)) - f(j)*dt*c/dx*2; % Inject source
   
    % DIPOLE SOURCE
    v(round(N/10)) =v(round(N/10)) - g(j)*dt/dx/rho*2; % Inject source
   
   % FINITE DIFFERENCE LEAP FROG. 4th order
   %     p(t+1/2,x) - p(t-1/2,x)              1/24u(t,x+3/2) - 9/*8u(t,x+1/2) + 9/8u(t,x-1/2) -1/24u(t,x-3/2)
   % (1) -----------------------    = -rhoc^2 ---------------------------------------------------------------
   %                dt                                                    dx
   %
   %     u(t+1,x+1/2) - u(t,x+1/2)            1/24p(t+1/2,x+2)- 9/8p(t+1/2,x+1)+ 9/8p(t+1/2,x) -1/24p(t+1/2,x-1)
   % (2) -------------------------  = -1/rho  ------------------------------------------------------------------
   %                dt                                                    dx
   % --------> Taylor series implementation as:
   p(XX+1) = -dt/dx*rho*c^2*( -1/24*v(XX-1) + 9/8*v(XX) - 9/8*v(XX+1) + 1/24*v(XX+2) ) + p(XX+1);
   v(XX)   = -dt/dx/rho*    ( -1/24*p(XX-1) + 9/8*p(XX) - 9/8*p(XX+1) + 1/24*p(XX+2) ) + v(XX  );
 
    % Optimal coefficient implementation as:
%     p(XX+1) = -dt/dx*rho*c^2*( co(2)*v(XX-1) + co(1)*v(XX) - co(1)*v(XX+1) -co(2)*v(XX+2) ) + p(XX+1);
%     v(XX)   = -dt/dx/rho*    ( co(2)*p(XX-1) + co(1)*p(XX) - co(1)*p(XX+1) -co(2)*p(XX+2) ) + v(XX  );
   
    %======================================================================
    % Record synthetic seismogram
    %======================================================================
    record(j) = p(round(3*N/5));
   
   
    %======================================================================
    % Plot
    %======================================================================
    plot(x,p),ylim([-1 1]),title(sprintf('Staggered grid, 2nd order time-Optimal space, time t=%0.2f',j*dt)),xlabel('offset [m]')
    hold on
    plot(x(round(3*N/5)),p(round(3*N/5)),'*')
    hold off
    drawnow
   

    %======================================================================
    % Store GIF
    %======================================================================
%     if mod(j,2)==0
%         fh=getframe(h);
%         im = frame2im(fh);
%         [imind,cm] = rgb2ind(im,256);
%           if gif == 1 
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
%           else 
%               imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%           end
%         gif = gif+1;
%     end
end

%% Recorded seismogram
figure(2)
plot(T,ITDT(record),'k')
title('recorded, with dispersion')
xlabel('time')