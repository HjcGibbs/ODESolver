function [tout,yout] = IRK4Solver(f, t, y0)
% INPUT: f(t,y) is an anonymous function that defines
% the right-hand side of the ODE ydot = f(t,y)
% t =[t0 t1 ... tfinal] is a vector of grid points
% with length N
% y0=[a b c] is a column vector that contain the
% initial values x(0)=a, y(0)=b, z(0)=c.
% OUTPUT:tout is a column vector of grid points.
% yout is an 3 x N matrix containing the solution
% at different grid points.

% Harry Gibbs       z3337442
% Sebastian Blefari z3416129
% Ben Madafiglio    z3460922

tout = t;
yout = [t;t;t];
yout(:,1) = y0;

% Define tableau (two-stage of order 4)
c = [1/2 - sqrt(3)/6, 1/2 + sqrt(3)/6];
b = [1/2, 1/2];
A = [1/4, 1/4 - sqrt(3)/6; 1/4 + sqrt(3)/6, 1/4];


for k = 1:(length(t) - 1)
    
    h = t(k+1)-t(k);
    yk = yout(:, k);

    XI1 = @(xi1, xi2) yk + h*A(1, 1)*f(t(k) + c(1)*h, xi1) + h*A(1, 2)*f(t(k) + c(2)*h, xi2) - xi1;
    XI2 = @(xi1, xi2) yk + h*A(2, 1)*f(t(k) + c(1)*h, xi1) + h*A(2, 2)*f(t(k) + c(2)*h, xi2) - xi2;
    
    p = @(p) [XI1(p(:,1),p(:,2)),XI2(p(:,1),p(:,2))];
    
    x0 = [yk, (yk + h*f(t(k), yk))];
    
    options=optimset('Display','off','TolFun',1e-32,'TolX',1e-32);
    psi = fsolve(p,x0,options); 
   
    yout(:, k + 1) = yk ...
        + h*b(1)*f(t(k) + c(1)*h, psi(:,1)) ...
        + h*b(2)*f(t(k) + c(2)*h, psi(:,2));
end

end