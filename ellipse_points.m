function [xpts,ypts] = ellipse_points(a,b,N)
% Returns the x,y points along the arc of an ellipse with equal arc length
% between each point with the first point at (0,b)
% Input: a - the major axis length
%        b -  minor axis length
%        N - the number of points on the circumference
% Output: vector of x and y points, and a figure demonstrating the points
% on the ellipse
% Source: most everything can be found at the link below...
% https://math.stackexchange.com/questions/172766/calculating-equidistant-points-around-an-ellipse-arc
% Reject 1 as a value
if N == 1
    fprintf('Please enter N value greater than 1.\n')
    return
end
% parameters:
m = 1-(b^2/a^2);
% circumference
%
fun = @(t) (sqrt(1-m.*(t.^2))./sqrt(1-(t.^2)));
E = integral(fun,0,1);
% this is the circumference of the ellipse
c = 4*a*E;
% First Arc length from x=0, y=b
arc_length = c/N;
% guess may not work for every N, it is an assumption that the arc length
% will not be much different than the arc length of a circle and its
% corresponding angle of 2*pi/N
    %phis = zeros(N,1);
    xpts = zeros(N,1);
    y_sign = zeros(N,1);
    for i=1:N
        phi_actual = myinverse(i*arc_length,i*2*pi/N,m,a);
        % phis(i) = phi_actual;
        xpts(i) = a*sin(phi_actual);
        % change the y sign for the quadrant...
        if (i*2*pi/N > pi/2 && i*2*pi/N < 3*pi/2)
            y_sign(i) = -1;
        else
            y_sign(i) = 1;
        end
    end
        ypts = sqrt(b^2-(b^2/a^2).*xpts.^2);
        ypts = ypts.*y_sign;
    % Ellipse values
%     x_ = -a:.1:a;
%     x_e = [x_,x_];
%     y_ = sqrt(b^2-(b^2/a^2).*x_.^2);
%     y_e = [y_,-y_];
    %figure
%     plot(x_e,y_e,'r*')
%     hold on
%     plot(xpts,ypts,'b*','LineWidth',10)
%     title(sprintf('Ellipse with %d equidistant points',N))
end
% functions
% returns the approximate amplitude for the given arc length
function x = myinverse(arc_length, guess_amp,m,a)
  amp_fun = @(theta) sqrt(1-m.*sin(theta).^2);
  E_phi = @(phi) a*integral(amp_fun,0,phi);
  fun = @(phi) E_phi(phi)-arc_length;  % parameterized function
  x = fzero(fun,guess_amp);     % guess is your initial pretty good guess for x
end