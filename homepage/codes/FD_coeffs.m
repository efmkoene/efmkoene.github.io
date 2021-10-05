function [c,varargout] = FD_coeffs( M, order, error )
% c = FD_coeffs( M, order, error ) computes M optimal FD coefficients to
% approximate an derivative of given order, with a given maximum wavenumber
% error, using a least squares algorithm.
%
% [c,max_wavenumber] = FD_coeffs( M, order, error ) returns the maximum 
% wavenumber within the desired error range.
%
%    order: The derivative order. 
%           1-> 1st order d/dx    . M = c_1, c_2, ..., c_M
%           2-> 2nd order d^2/dx^2. M = c_0, c_1, ..., C_{M-1}
%    error: Typically 1e-4. 
%
% Note that for the 1st order, one wants to flip the sign for x_{i-M}
% points.
% 
%
% Based on Liu (2013): http://library.seg.org/doi/pdf/10.1190/geo2012-0480.1
% and Liu (2014): http://dx.doi.org/10.1093/gji/ggu032
% Though (for some reason?) the coefficients are very slightly different.
% 
% Erik Koene, 24-2-2017

% %%%%%%%%%% The optimization functions
if order == 1
    % The computed wavenumber
    cfn = @(m,n,b) sum( 4 * ( sin((m+1/2)*b)-2*(m+1/2)*sin(b/2) ) ...
                         .* ( sin((n+1/2)*b)-2*(n+1/2)*sin(b/2) ) ...
                                                                  )./length(b);
    
    % The desired wavenumber
    dfn = @(n,b) sum(   2 * ( sin((n+1/2)*b)-2*(n+1/2)*sin(b/2) ) ...
                         .* (             b -2*        sin(b/2) ) ...
                                                                  )./length(b);
    
    % The cost function
    cf  = @(b,x,M) max( x'* 2*sin(([2:M]-1/2)'*b) - b + (1-([3:2:(2*M-1)]*x))*2*sin(b/2) );

elseif order == 2
    % The computed squared wavenumber
    cfn = @(m,n,b) sum( 16*sin( m*b/2 ).^2 .* sin( n*b/2 ).^2 )./length(b);
    
    % The desired squared wavenumber
    dfn = @(n,b) sum( 4 * (b).^2 .* sin( n*b/2 ).^2 )./length(b);
    
    % The cost function
    cf  = @(b,x,M) max( x' * 4* sin( [1:M-1]' * b/2 ).^2 - b.^2 );
end

% %%%%%%%%%% Testing the error for different wavenumber ranges
db = 0.005;
wavenumbers = db:db:pi;
wavenumerror = [];
C = zeros(M-1,M-1);
d = zeros(M-1, 1 );
for bmax=wavenumbers
    % Wavenumber range
    b = [0:db:bmax];

    % System of equations for the tested wavenumber range
    for m=1:M-1
        for n=1:M-1
            C(n,m) = cfn(m,n,b);
        end
        d(m) = dfn(m,b);
    end

    % Least squares solution
    x=C\d;

    % Compute the maximum error
    wavenumerror = [wavenumerror ; cf(b,x,M)];
end

% %%%%%%%%%% Find the maximum wavenumber within the error bound
tmp = wavenumbers( abs(wavenumerror)<error );

% %%%%%%%%%% System of equations for the desired wavenumber range
b = [0:db:tmp(end)];
for m=1:M-1
    for n=1:M-1
        C(n,m) = cfn(m,n,b);
    end
    d(m) = dfn(m,b);
end

% %%%%%%%%%% Find least squares solution for desired wavenumber range
x=C\d;

% %%%%%%%%%% Generate list of coefficients
if order==1
    c(1)   = 1-([3:2:2*M-1])*x;
elseif order==2
    c(1)   = -2*sum(x);
end
c(2:M) = x;

for k = 1:max(nargout,1)-1
    varargout{k} = tmp(end);
end

end