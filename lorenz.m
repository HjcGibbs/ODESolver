function y = lorenz(~,x,sigma,rho,beta)
    % INPUT: t is a a real value indicating time
    % x is a column vector of size 3 x1
    %        sigma, rho, beta are parameters of the Lorenz
    %        equations
    % OUPUT: y is a column vector of size 3 x 1 that gives
    %        the right hand side of the Lorenz equations
    
    % Harry Gibbs       z3337442
    % Sebastian Blefari z3416129
    % Ben Madafiglio    z3460922
    
    a = x(1,1);
    b = x(2,1);
    c = x(3,1);
    y = [sigma*(b - a);a*(rho - c) - b; a*b - beta*c];
end
    
    
    