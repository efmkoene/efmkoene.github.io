close all; clear all; clc;

% Process gif with Bash command: gifsicle -O3 -d3 -b --lossy=80 --colors 40 --crop 33,20+-25x-20 time_dispersion.gif

%% Model setup
c = 2000;                  % [m/s] wave velocity
dx= 5;                     % [m]   grid spacing
dt= .0015;                 % [s]   time sampling
r = c*dt/dx                % [-]   CFL (Courant number)
x = 0:dx:2000;             % [m]   cell locations
N = length(x);             % [-]   number of cells
T_l = 0:dt:1;              % [s]   time vector long step
short = 10;                % [-]   Subsampling of dt (1/short * dt)
T_s = 0:dt/short:1;        % [s]   time vector short step
M = length(T_l);           % [-]   number of time steps
xs = round(N/10);          % [-]   Source injection location
xr = round(3*N/5);         % [-]   Receiver seismogram location

% Ricker wavelet derivative source function
ricker_d = @(fc,t) 2*pi^2*fc^2*t .* (2*pi^2*fc^2*t.^2 - 3) .* exp( -pi^2*fc^2.*t.^2 );
fc = 35;
d  = 0.04;
fs_l = ricker_d(fc,T_l-d);
fs_s = ricker_d(fc,T_s-d);

% Initialize the grid
[u1_l, u2_l, uj_l, u1_s, u2_s, uj_s ] = deal( zeros(N,1) );
% u1 @ t-dt
% u2 @ t
% uj @ t+dt

% Seismogram init
record = zeros(M,2);

% GIF INIT
h=figure(1); set(gcf,'Position',[2700 200 562 557]); axis tight;
filename='../images/time_dispersion.gif';
gif = 1;

X = 2:N-1; % Update cells
for j=2:M
    %======================================================================
    % Long dt
    %======================================================================
    % Source injection in update
    u2_l(xs) =u2_l(xs) + fs_l(j)*dt^2/dx*c*2; % Inject source
    
    % Pseudospectral medium update ( -k^2*u in frequency domain )
    k = ([(pi/dx/N)*(0:N-1) -pi/dx+(pi/dx/N)*(0:N-1)]);
    u = fft(u2_l(:),2*N);
    tmp = ifft(-k.^2' .* u,'symmetric');
    
    % Update step
    uj_l(X) = dt^2*c^2*tmp(X) + 2*u2_l(X) - u1_l(X);
    
    % Prepare for next update
    u1_l = u2_l;
    u2_l = uj_l;

    %======================================================================
    % Short dt
    %======================================================================
    for l=2:short+1
        % Source injection in update
        u2_s(xs) =u2_s(xs) + fs_s((j-2)*short+l)*dt/M*N*2/short^2; % Inject source

        % Pseudospectral medium update ( -k^2*u in frequency domain )
        k   = ([(pi/dx/N)*(0:N-1) -pi/dx+(pi/dx/N)*(0:N-1)]);
        u  = fft(u2_s(:),2*N);
        tmp = ifft(-k.^2' .* u,'symmetric');

        % Update step
        uj_s(X) = (dt/short)^2*c^2*tmp(X) + 2*u2_s(X) - u1_s(X);

        % Prepare for next update
        u1_s = u2_s;
        u2_s = uj_s;
    end
    
    
    %======================================================================
    % Record synthetic seismogram
    %======================================================================
    record(j,:) = [uj_l(xr), uj_s(xr)];
   
    %======================================================================
    % Plot
    %======================================================================
    subplot(3,1,1)
    plot(x,uj_l),ylim([-1 1]),title(sprintf('Long dt=%0.2f ms; time t=%0.2f s',dt*1e3,j*dt)),xlabel('offset [m]')
    hold on
    plot(x(xr),uj_l(xr),'*')
    hold off
    
    subplot(3,1,2)
    plot(x,uj_s),ylim([-1 1]),title(sprintf('Short dt=%0.2f ms; time t=%0.2f s',dt/short*1e3,j*dt)),xlabel('offset [m]')
    hold on
    plot(x(xr),uj_s(xr),'*')
    hold off
    
    subplot(3,1,3)
    plot( T_l, record ), title(sprintf('Seismogram at time t=%0.2f s',j*dt))
    hold on
    plot( T_l(j), uj_s(xr),'*' )
    plot( T_l(j), uj_l(xr),'*' )
    hold off
    ylim([-1 1])
    xlabel('Time [s]')
    legend('Long \Deltat','Short \Deltat','Location','NorthWest')
    
    drawnow
    
    %======================================================================
    % Store GIF
    %======================================================================
    if mod(j,2)==0
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
plot( T_l, record )
legend('Long \Deltat','Short \Deltat')
title('Recorded synthetic seismogram')