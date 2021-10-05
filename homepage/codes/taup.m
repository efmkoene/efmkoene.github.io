function [out] = taup( input_data, dt, offset, p, epsilon, mode, cut)
% Erik Koene, April 2016.
%
% Inputs:
%   input_data [MAT] = traces(nt, nx) or radon transform(nt, np)
%   dt         [SCA] = sample spacing dt
%   offset     [VEC] = trace offsets in x wrt. near-offset, i.e the intercept point relative to slowness lines 
%   ps         [VEC] = slowness vector for example linspace(-1/vp, 1/vp, nx)
%   epsilon    [SCA] = stabilisation factor (~1E-3 works well)
%   mode       [SCA] = -1: inverse, 1: forward
%   cut        [SCA] = nt, the number of samples in one trace (i.e length(P,1) ) 
%
% Example usage:
%   P_tp  = taup( P   , dt, xs, linspace(-1/vp, 1/vp, 400), 1E-3,  1,    0      )
%   P_rec = taup( P_tp, dt, xs, linspace(-1/vp, 1/vp, 400),  0  , -1, size(P,1) )
%   fprintf(1,'Amplitudes P  correct up to %f %%\n', 100*max(max(abs(P-P_rec  )))/max(max(abs(P ))));

    % Input dimensions
    nt = size(input_data,1);
    np = length(p);
    nx = length(offset);
    
    % Force shape of offset and p vectors to be orthogonal
    p      = reshape(p     , np, 1 );
    offset = reshape(offset, 1 , nx);
    px     = p * offset;
    
    % Fourier transform input data over TIME
    if mode==1            
        nfft   = 2^( nextpow2(nt)+1 );
    else
        nfft   = nt;
    end
    D_fft  = fft( input_data, nfft, 1 );
    
    % Extend to zero-padded data domain
    tmax = dt * nfft;
    df   = 1/tmax;
    
    % Initialize zeros in fft'd data
    if mode==1
        out_fft = zeros( nfft, np );
    elseif mode==-1
        out_fft  = zeros( nfft, nx );
    end
       
    % Loop over frequencies, for now take entire frequency range 0:f_nyq
    % (other programs have options to set f_min and f_max)
    ifmin = 1;
    ifmax = nfft/2 + 1;
    for ifx = ifmin:ifmax
       w = 2*pi*(ifx-1)*df; % omega
       
       % Create matrix L x->p (takes longest time. speed up???)
       L = exp( 1i * w * px );
       
       % Least squares forward transform
       if mode==1
           if np <= nx
               LLH = L*L';
               LLH = LLH + epsilon*max(abs(diag(LLH)))*eye(np); 
               out_fft(ifx,:) = LLH \ (L * D_fft(ifx,:).') ;  % Invert xp = (L*L')^-1 * L * xin
           elseif np > nx
               LHL = L'*L;
               LHL = LHL + epsilon*max(abs(diag(LHL)))*eye(nx);
               out_fft(ifx,:) = L * (LHL \ D_fft(ifx,:).');   % Invert xp = L * (L'*L)^-1 * xin
           end
       % Fast inverse transform
       elseif mode==-1
          out_fft(ifx,:) = L' * D_fft(ifx,:).' ;
       end
    end
    
    % Transform back to time-domain
    if mode==1
        out_fft(1,:) = 0.5*out_fft(1,:);
        
        out = real( ifft( out_fft , [], 1) );
        out = 2*out;
    elseif mode==-1
        out_fft(1,:) = 0.5*out_fft(1,:);
        
        out  = real( ifft( out_fft, [], 1) );
        out  = 2*out(1:cut, :);
    end
end