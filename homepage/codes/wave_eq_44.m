close all; clear all; clc;

% Process gif with Bash command: gifsicle -O3 -d3 -b --lossy=80 --colors 40 --crop 33,20+-25x-20 space_dispersion.gif

%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
dt= 0.9*(dx/c);            % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
x = 0:dx:2000;             % [m]   cell locations
N = length(x);             % [-]   number of cells
T = 0:dt:1;                % [s]   time vector long step
M = length(T);             % [-]   number of time steps
xs = round(N/10);          % [-]   Source injection location
xr = round(3*N/5);         % [-]   Receiver seismogram location

% Ricker wavelet derivative source function
ricker_d = @(fc,t) 2*pi^2*fc^2*t .* (2*pi^2*fc^2*t.^2 - 3) .* exp( -pi^2*fc^2.*t.^2 );
fc = 35;
d  = 0.04;
fs = ricker_d(fc,T-d);

% Initialize the grid
[ u1, u2, uj ] = deal( zeros(N,1) );
% u1 @ t-dt
% u2 @ t
% uj @ t+dt

% Seismogram init
record = zeros(M,2);

% GIF INIT
h=figure(1); set(gcf,'Position',[2700 200 562 557]); axis tight;
filename='../images/space_dispersion.gif';
gif = 1;

X = 3:N-2; % Update cells
X5= 5:N-4;

for j=2:M
    %======================================================================
    % Simulation with Lax-Wendroff correction
    %======================================================================
    % Source injection in update
    u2(xs) =u2(xs) + fs(j)*dt/M*N*2; % Inject source
        
    % Update step
    du = ( -1/12*u2(X-2) + 4/3*u2(X-1) - 5/2*u2(X) ...
           -1/12*u2(X+2) + 4/3*u2(X+1) );
    du2= ( -1/12*du(1:end-4) + 4/3*du(2:end-3) - 5/2*du(3:end-2) ...
           -1/12*du(4:end-1) + 4/3*du(5:end) );
    uj(X) = r^2*du + 2*u2(X) - u1(X);     % Normal update
    uj(X5)= r^4/16*du2 + uj(X5);          % LW correction
    
    % Prepare for next update
    u1 = u2;
    u2 = uj;
    
    %======================================================================
    % Record synthetic seismogram
    %======================================================================
    record(j,:) = [uj(xr)];
   
    %======================================================================
    % Plot
    %======================================================================
    subplot(2,1,1)
    plot(x,uj),ylim([-1 1]),title(sprintf('Short dx=%0.2f m; time t=%0.2f s',dx,j*dt)),xlabel('offset [m]')
    hold on
    plot(x(xr),uj(xr),'*')
    hold off
    
    subplot(2,1,2)
    plot( T, record ), title(sprintf('Seismogram at time t=%0.2f s',j*dt))
    hold on
    plot( T(j), uj(xr),'*' )
    hold off
    ylim([-1 1])
    xlabel('Time [s]')
    legend('Short \Deltat','Location','NorthWest')
    
    drawnow
    
    %======================================================================
    % Store GIF
    %======================================================================
%     if mod(j,6)==0
%         f=getframe(h);
%         im = frame2im(f);
%         [imind,cm] = rgb2ind(im,256);
%           if gif == 1 
%               imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
%           else 
%               imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%           end
%         gif = gif+1;
%     end

end

%% Plot
plot( T, record )
legend('Long \Deltax','Short \Deltax')
title('Recorded synthetic seismogram')
ylim([-1 1])