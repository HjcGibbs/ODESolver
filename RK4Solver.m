function [tout,yout] = RK4Solver(f,t,y0);
    % INPUT: f(t,y) is an anonymous function that defines
    %        the right-hand side of the ODE ydot = f(t,y)
    %        t =[t0 t1 ... tfinal] is a vector of grid points
    %           with length N
    %        y0=[a b c] is a column vector that contain the
    %        initial values x(0)=a, y(0)=b, z(0)=c.
    % OUTPUT:tout is a column vector of grid points.
    %        yout is an 3 x N matrix containing the solution
    %        at different grid points.
    
% Harry Gibbs       z3337442
% Sebastian Blefari z3416129
% Ben Madafiglio    z3460922
    
    tout = t;
    yout = [t;t;t];
    yout(:,1) = y0;
    
 for k = 1:length(t)-1
     
     h = t(k+1)-t(k);
     
     currY = yout(:, k);
     currT = t(k);
     
     k1 = f(currT, currY);
     k2 = f(currT + h/2, currY + h*k1/2);
     k3 = f(currT + h/2, currY + h*k2/2);
     k4 = f(currT + h, currY + h*k3);
     yout(:, k+1) = currY + h*(k1 + 2*k2 + 2*k3 + k4)/6;
 end
    
 end