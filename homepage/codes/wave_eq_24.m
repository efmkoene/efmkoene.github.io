close all; clear all; clc;

% Process gif with Bash command: gifsicle -O3 -d3 -b --lossy=80 --colors 40 --crop 33,20+-25x-20 space_dispersion_24.gif

%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
dt= .0005;                 % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
x = 0:dx:2000;             % [m]   cell locations
short = 3;                 % [-]   Subsampling of dt (1/short * dt)
x2= 0:dx/short:2000;       % [m]   shortened cell locations
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
[ u1_l, u2_l, uj_l ] = deal( zeros(N,1) );
[ u1_s, u2_s, uj_s ] = deal( zeros(short*(N-1)+1,1) );
% u1 @ t-dt
% u2 @ t
% uj @ t+dt

% Seismogram init
record = zeros(M,2);

% GIF INIT
h=figure(1); set(gcf,'Position',[2700 200 562 557]); axis tight;
filename='../images/space_dispersion_24.gif';
gif = 1;

X = 3:N-2; % Update cells
X2= 3:short*(N-1)-1; % For short

for j=2:M
    %======================================================================
    % Large dx
    %======================================================================
    % Source injection in update
    u2_l(xs) =u2_l(xs) + fs(j)*dt/M*N*2; % Inject source
        
    % Update step
    uj_l(X) = r^2 * (-1/12*u2_l(X-2) + 4/3*u2_l(X-1) - 5/2*u2_l(X) ...
                     -1/12*u2_l(X+2) + 4/3*u2_l(X+1))+ 2*u2_l(X) - u1_l(X);
    
    % Prepare for next update
    u1_l = u2_l;
    u2_l = uj_l;

    %======================================================================
    % Short dx
    %======================================================================
    % Source injection in update
    u2_s((xs-1)*short+1) =u2_s((xs-1)*short+1) + fs(j)*dt/M*N*short*2; % Inject source
        
    % Update step
    uj_s(X2) = r^2*short^2 * (-1/12*u2_s(X2-2) + 4/3*u2_s(X2-1) - 5/2*u2_s(X2) ...
                              -1/12*u2_s(X2+2) + 4/3*u2_s(X2+1))  + 2*u2_s(X2) - u1_s(X2);
    
    % Prepare for next update
    u1_s = u2_s;
    u2_s = uj_s;
    
    %======================================================================
    % Record synthetic seismogram
    %======================================================================
    record(j,:) = [uj_l(xr), uj_s(xr*short)];
   
    %======================================================================
    % Plot
    %======================================================================
    subplot(3,1,1)
    plot(x,uj_l),ylim([-1 1]),title(sprintf('4th order in space, long dx=%0.2f m; time t=%0.2f s',dx,j*dt)),xlabel('offset [m]')
    hold on
    plot(x(xr),uj_l(xr),'*')
    hold off
    
    subplot(3,1,2)
    plot(x2,uj_s),ylim([-1 1]),title(sprintf('4th order in space, short dx=%0.2f m; time t=%0.2f s',dx/short,j*dt)),xlabel('offset [m]')
    hold on
    plot(x(xr),uj_s(xr*short),'*')
    hold off
    
    subplot(3,1,3)
    plot( T, record ), title(sprintf('Seismogram at time t=%0.2f s',j*dt))
    hold on
    plot( T(j), uj_s(xr*short),'*' )
    plot( T(j), uj_l(xr),'*' )
    hold off
    ylim([-1 1])
    xlabel('Time [s]')
    legend('Long \Deltax','Short \Deltax','Location','NorthWest')
    
    drawnow
    
    %======================================================================
    % Store GIF
    %======================================================================
    if mod(j,6)==0
        f=getframe(h);
        im = frame2im(f);
        [imind,cm] = rgb2ind(im,256);
          if gif == 1 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end
        gif = gif+1;
    end

end

%% Plot
plot( T, record )
legend('Long \Deltax','Short \Deltax')
title('Recorded synthetic seismogram')