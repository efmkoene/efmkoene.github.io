%% n-th order space / 2nd order in time FD example
% %%%%%%%%%%%%%%%%% Design the FD operator along with DT and DX
% INPUTS
fc   = 25;    % [Hz]  central wavelet frequency
minc = 2000;  % [m/s] minimum velocity
maxc = 3000;  % [m/s] maximum velocity
order= 10;     % FD operator order (2*order-1 point stencil)

% FD OPERATOR DESIGN
error   = 5e-5; % Fixed phase error level (1e-4 or 1e-5 are common).
[fdc,e] = FD_coeffs( order, 1, error );  % FD coefficients
s = sum( abs(fdc) )^-1;                  % (CFL) Stability number

% LIU
% fdc = [1.264748  -0.1331606   0.04296909  -0.0181897   0.008861071  -0.004347073    0.002076101  -0.0009164925   0.0003437446  -0.00007874250]


% ESTABLISH DT and DX
dx = minc/(fc*3)/(2*pi/e); % dx -> Pts./wavelength, based on max-Ricker-wavelet-freq: 2.5*fc Hz
dt = s*dx/maxc;              % dt -> Chosen by CFL stability number
fprintf('dx=%f,  dt=%f\n',dx,dt)

% %%%%%%%%%%%%%%%%% Init the simulation itself
% SIMULATION DOMAIN
x = 0:dx:5000;             % [m]   cell locations
N = length(x);             % [-]   number of cells
et= 1.5;                  % [s]   end time
T = 0:dt:et;               % [s]   time vector
M = length(T);             % [-]   number of time steps
% Initialize the grid
p1 = zeros(N,1); % 'u 2 old'
v2 = zeros(N,1); % 'u new'
% Initialize the seismogram
record = zeros(M,1);

% RICKER WAVELET SOURCE
t0 = 1.5/fc;               % [s]   delay wavelet
tau=pi*(T-t0)*fc;          % [-]   Ricker wavelet argument
fs = (1 - 2*tau.^2) .* exp( -tau.^2 );       % Ricker wavelet

% ADD DISPERSION TO SOURCE
fs = FTDT(fs);

% VARIABLE VELOCITY MODEL
c = ones(N,1)*minc;
c(end/2:end)=maxc;
r = c*dt/dx;               % [-]   CFL (Courant number)
fprintf('Max CFL number %f <?= %f\n',max(r),s)
fprintf('Pts/wavelength %f >?= %f)\n',min(c)/(fc*2.5)/dx, 2*pi/e)

% %%%%%%%%%%%%%%%%% The simulation
% CHEAP LEFT BOUNDARY CANCELATION
cancel = fs( abs(fs) > (max(fs)*0.001));
cancel = max(find(fs == cancel(end) ))+1;   % Last injection index
% cancel = M;

% SOURCE INJECTION POINT
xinj = round(order+(cancel*dt*minc)/(2*dx));

% POINT UPDATE RANGE
X = order+1:N-order;

% SIMULATION LOOP
for j=2:M
    % SOURCE INJECTION
    p1(xinj) =p1(xinj) +  fs(j)*sqrt(dt/M)*N*2;
    
    % TIME INTEGRATION (TIME STEPPING)
    % First component
    for i=1:order
       p1(X+1) = p1(X+1) + c(X).^2*dt/dx .*fdc(i).*( v2(X+1-i) - v2(X+i) );
    end
    % Second component
    for i=1:order
       v2(X  ) = v2(X  ) +         dt/dx .*fdc(i).*( p1(X+1-i) - p1(X+i) );
    end
    
    % CANCEL THE LEFT-GOING WAVE
    if j > cancel
        p1(1:xinj)=0;
        v2(1:xinj)=0;
    end
    
    % PLOT THE RESULTS
    plot(x,p1),ylim([-1 1]),title(sprintf('time t=%0.2f',j*dt)),xlabel('offset [m]')
    hold on
    plot(x(round(3*N/5)),p1(round(3*N/5)),'*')
    hold off
    drawnow
    
    % RECORD THE SEISMOGRAM
    record(j) = p1(round(3*N/5));
end

% SAVE THE SEISMOGRAM
coarse_seismogram = record;

figure(2)
plot(T,coarse_seismogram)
title('Recorded trace')
xlabel('Time [s]')
legend('Coarse in time')
% xlim([0.42 0.59])
