function [ itdt ] = ITDT( f, varargin )
%ITDT Inverse Time Dispersion Transform
%   This function artificially REMOVES time dispersion from a vector of samples
%   in time. It models the phase shift error. It uses the fact that every
%   sample equals 1 time step (assuming equidistant sampling), removing the dt dependency.
%   f        [VEC]        measurements
%   varargin [INT] (opt.) subsampling rate
%  USE: 
%      ITDT( data    ); % Removes dispersion from entire trace, assuming 1 sample equals 1 simulation timestep.
%      ITDT( data, 4 ); % Removes dispersion from entire trace, assuming 1 sample equals 4 simulation timesteps.
%
% 10-2016 EK: Basic version (Wang & Xu, 2015)
% 1 -2017 EK: Various optimisations.

% Amount of subsampling
if length(varargin)==1
    ss = varargin{1};
else
    ss = 1;
end

% Zero-padd to twice the size
f = [        f(:)       ; 
    zeros(length(f), 1 )];

% The phase shift function
fn = @(omega) 2*asin( omega/2 );

% FFT
nf = length(f);

% Take altered Fourier Transform
IFTDTuf= exp( -1i * fn([0:floor(nf/pi)]*2*pi/nf/ss)'*ss * [0:nf-1]) * f(:); % Apply ITDT in an altered DFT
IFTDTuf(ceil(nf/pi):nf)= 0; % Block out omega*dt>2 ==> f>1/pi

% Take back to the time domain
itdt   = ifft( IFTDTuf, 'symmetric');                         % Bring back to the time domain 
itdt   = itdt(1:end/2);                                       % Remove zero-padded entries 

end