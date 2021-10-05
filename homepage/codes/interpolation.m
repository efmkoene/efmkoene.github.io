% Function
y = @(x) 1./(1+x.^2);

% Fine spaced solution
x0 = linspace(-5,5,200);
y0 = y(x0);

%==========================================================================
% Standard polynomial interpolation
%==========================================================================
figure(1); 
plot( x0, y0, 'k', 'LineWidth',2 ); % True solution
ax=gca;
j=1;
for N=[4 7 16]
    x=linspace(-5,5,N);      % Subsample at N points
    p = polyfit(x,y(x),N-1); % Create polynomial of degree N-1
    
    hold on;
        
        ax.ColorOrderIndex = j+2;
        plot( x0, polyval(p,x0), '-.','Linewidth',1.5 );
        ax.ColorOrderIndex = j+2;
        h=plot( x , y(x) , 'o');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    hold off % Plot
    j=j+1;
end
grid minor

%==========================================================================
% Plot
%==========================================================================
ylim([-.5 1])
title('Approximating a function with regular nodes')
legend('Function: 1/(1+x^2)','3rd degree polynomial','6th degree polynomial','15th degree polynomial','Location','South')
set(gcf,'Position',[942 618 575 280])


%==========================================================================
% Chebyshev (cosine) interpolation
%==========================================================================
figure(2);
plot( x0, y0, 'k', 'LineWidth',2 ); % True solution
ax=gca;
ax.ColorOrderIndex=1;
j=3;
for N=[4 7 16]
    x=linspace(-5,5,N);      % Subsample at N points
    x=5*cos((x+5)/10*pi);    % Get Gauss-Lobatto points
    p = polyfit(x,y(x),N-1); 
    
    hold on;
        ax.ColorOrderIndex = j;
        plot( x0, polyval(p,x0), '-.','Linewidth',1.5 );
        ax.ColorOrderIndex = j;
        h=plot( x , y(x) , 'o');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    hold off % Plot
    j=j+1;
end
grid minor
ylim([-.5 1])
title('Approximating a function with Chebyshev nodes')
legend('Function: 1/(1+x^2)','3rd degree polynomial','6th degree polynomial','15th degree polynomial','Location','South')
set(gcf,'Position',[942 254 575 280])