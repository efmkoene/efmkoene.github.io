% Cosine series approximation
fo = @(x) 1 - x.^2/factorial(2);                                           % 2nd order
so = @(x) 1 - x.^2/factorial(2) + x.^4/factorial(4);                       % 4th order
to = @(x) 1 - x.^2/factorial(2) + x.^4/factorial(4) - x.^6/factorial(6);   % 6th order

x=[0:0.01:2*pi];                                                           % x-range for 1 cosine

% Plot
figure(1)
plot( x, cos(x), x, fo(x), x, so(x), x, to(x) )
ylim([-1 1])
title('Cosine series approximation')
legend('cos(x)','2nd order Taylor expansion','4th order Taylor expansion','6th order Taylor expansion')

set(gcf,'Position',[680 823 680 280])