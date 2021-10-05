%% Model setup
%==========================================================================
% MODEL SETUP
%==========================================================================
c = 2000;                 % [m/s] velocity
dx= 5;                    % [m]   grid spacing x
dy= dx;                   % [m]   grid spacing y
r = 0.6;                  % [-]   Courant (CFL) number

x = 0:dx:2000;            % [m]   cell locations x
y = 0:dy:2000;            % [m]   cell locations y

%==========================================================================
% AUTOMATIC FTDT SETTINGS
%==========================================================================
N = length(x);            % [-]   number of cells

dt= r*dx/c;               % [s]   Time step sampling based on CFL
T = 0:dt:0.95;            % [s]   time vector
M = length(T);            % [-]   number of time steps

%==========================================================================
% SET SOURCE/RECEIVER LOCATIONS
%==========================================================================
is = ceil(N/2);
isx= round(x(end)/4/dx);
xr = round(N/2);
yr = round(N/2);

%==========================================================================
% SOURCE WAVELET DEFINITION
%==========================================================================
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
fs = ricker(30,T-0.06);

%==========================================================================
% 45 deg. filter (using Hilbert transform)
%==========================================================================
% fs_h  = hilbert(fs);
% phase = -pi/4;
% fs    = cos(phase)*real(fs_h) - sin(phase)*imag(fs_h);

%==========================================================================
% Half-order derivative (multiply with sqrt(-i*omega) in freq. domain)
%==========================================================================
% fs_f  = fft(fs);
% omega = 2*pi* [1:round(M/2)-1]/2/dt/M;         % omega=2*pi*f
% fs_f(2:round(M/2)) = fs_f(2:round(M/2)) .* sqrt(-i*omega);   % Applying the filter
% fs_f(M/2+1:end) = 0;
% fs = ifft(fs_f,'symmetric');


%==========================================================================
% Fix Bessel function behaviour: G_2d = -1/4 i H_0^2( omega*t ).
%==========================================================================
fs_f  = fft(fs);
omega = 2*pi* [1:round(M/2)-1]/2/dt/M;         % omega=2*pi*f
fs_f(2:round(M/2)) = fs_f(2:round(M/2)) .* besselk(0,i*omega);   % Applying the filter
fs_f(M/2+1:end) = 0;
fs = ifftshift(ifft(fs_f,'symmetric'));


%==========================================================================
% SEISMOGRAM SETUP
%==========================================================================
record = zeros(M,1);

%==========================================================================
% INIT FIELDS
%==========================================================================
[px,py,vx,vy] = deal( zeros( N, N ) );

%==========================================================================
% COMPLICATE VELOCITY MODEL
%==========================================================================
c  = ones(N, N) * c;
c(:, 1:round(N/4)) = 1500;

%==========================================================================
% UPDATE REGION
%==========================================================================
XX= 2:N-2;

%==========================================================================
% RUN SIMULATION
%==========================================================================
figure(1)
for j=2:M
    % Source injection
    px( isx, is ) = px(isx,is) + 2*pi*dt/dx^2*c(isx,is)^2*fs(j);
    py( isx, is ) = py(isx,is) + 2*pi*dt/dx^2*c(isx,is)^2*fs(j);
                           
    % FINITE DIFFERENCE LEAP FROG. 4th order
    %     p(t+1/2,x) - p(t-1/2,x)          u(t,x+1/2) - u(t,x-1/2)
    % (1) -----------------------    = c^2 -----------------------
    %                dt                               dx
    %
    %     u(t+1,x+1/2) - u(t,x+1/2)        p(t+1/2,x+1) - p(t+1/2,x)
    % (2) -------------------------  =     -------------------------
    %                  dt                               dx
    px(XX+1,XX+1) = dt/dx*c(XX+1,XX+1).^2.*(  -1/24*vx(XX-1,XX) + 9/8*vx(XX,XX) - 9/8*vx(XX+1, XX  ) + 1/24*vx(XX+2,XX) ) + px(XX+1,XX+1);
    py(XX+1,XX+1) = dt/dy*c(XX+1,XX+1).^2.*(  -1/24*vy(XX,XX-1) + 9/8*vy(XX,XX) - 9/8*vy(XX, XX+1  ) + 1/24*vy(XX,XX+2) ) + py(XX+1,XX+1);   
   
    vx(XX,XX)     = dt/dx*    ( -1/24*px(XX-1,XX+1) + 9/8*px(XX  ,XX+1) - 9/8*px(XX+1,XX+1) + 1/24*px(XX+2,XX+1) ...
                              + -1/24*py(XX-1,XX+1) + 9/8*py(XX  ,XX+1) - 9/8*py(XX+1,XX+1) + 1/24*py(XX+2,XX+1) ) + vx(XX,XX);
    vy(XX,XX)     = dt/dy*    ( -1/24*px(XX+1,XX-1) + 9/8*px(XX+1,XX  ) - 9/8*px(XX+1,XX+1) + 1/24*px(XX+1,XX+2) ...
                              + -1/24*py(XX+1,XX-1) + 9/8*py(XX+1,XX  ) - 9/8*py(XX+1,XX+1) + 1/24*py(XX+1,XX+2) ) + vy(XX,XX);
   
    % Record seismogram
    record(j) = px(xr,yr)+py(xr,yr);
    
    % Plot                      
    subplot(3,1,[1 2])
    imagesc(x,y,px+py),title(sprintf('time t=%0.2f',j*dt))
    caxis([-1 1])
    hold on
    plot(x(xr),y(yr),'r*')
    hold off
    axis equal
    subplot(3,1,3)
    plot(T,record)
    hold on
    plot( T(j), record(j), '*')
    hold off
    drawnow
end

%% Plot seismogram
figure(2)
plot(T,record)
title('Recorded seismogram')