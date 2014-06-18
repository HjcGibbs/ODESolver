function [ tout, yout ] = EulerSolver( t, f, y0 )
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
    yout(:,k+1) = h*f(t(k),yout(:,k)) + yout(:,k);
end

end