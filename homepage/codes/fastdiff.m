function Dx=fastdiff(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fastdiff.m
% 
% Numerically differentiate a function which has been evaluated at the 
% Chebyshev-Gauss-Lobatto points and stored in the vector x.
%
% Equivant to D*x where D is the Chebyshev collocation matrix, but much 
% faster when x contains many elements. The differentiation is performed
% in the spectral domain, by inverting the integral operator on a subspace.
% Since the integral operator has tridiagonal form, the entire
% differentiation process is O(N*log(N)).
%
% Reference on Chebyshev integration operators:
% 
% "Integration preconditioners for differential operators in spectral tau-
% methods" by E. A. Coutsias, T. Hagstrom, J. S. Hesthaven, and D. Torres. 
% Houston Journal of Mathematics p.21-38 (1996)
% Available here: http://www.math.unm.edu/~vageli/papers/ico.ps
%
% Written by: Greg von Winckel - 09/14/05
% Contact: gregvw(at)math(dot)unm(dot)edu
% URL: http://math.unm.edu/~gregvw
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N1=length(x); N=N1-1; k=1:N1; N2=N1+1;
band=[0 0 1./(2*(1:N2))]';
B=spdiags([band(k+2) -band(k)],[-1,1],N1,N1); 
B(2,1)=1; 
q=real(ifft([x(k);x(N:-1:2)]));
y=([q(1); 2*q(2:N); q(N1)]);
Dy=[B(2:N1,1:N)\y(2:N1);0];
q=fft([Dy(1); [Dy(2:N)/2; Dy(N1); Dy(N:-1:2)/2]]); 
Dx=real(q(k));