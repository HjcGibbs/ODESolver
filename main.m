% Harry Gibbs       z3337442
% Sebastian Blefari z3416129
% Ben Madafiglio    z3460922

sigma = 10; rho = 15; beta = 8/3;
f = @(t,x) lorenz(t,x,sigma,rho,beta);
k = 5;   % or 2,3,4,5,6
h = 10^(-k);
tfinal = 10;
t = [0:h:tfinal];
y0 = [-1;3;4];

solve = input('What method do you want to use? Euler, RK4 or IRK4. ');
if strncmpi(solve,'Euler',5)
    [tout, Y] = EulerSolver(t,f,y0);
elseif strncmpi(solve, 'rk4',3)
    [tout , Y] = RK4Solver(f,t,y0);
elseif strncmpi(solve,'irk4',4)
    [tout , Y] = IRK4Solver(f,t,y0);
else
    fprintf('You did not select a method available, try again.\n')
    return
end

options = odeset('RelTol',3.1e-14,'AbsTol',1e-16);
[tmout , Ym] = ode45(f,t,y0,options);


subplot(2,2,1)
plot(tmout,Ym)
title('ode45')
legend('x', 'y', 'z');
xlabel('t')
subplot(2,2,2)
plot3(Ym(:,3), Ym(:,2), Ym(:,1)) 
xlabel('z')
ylabel('y')
zlabel('x')
title('ode45')

subplot(2,2,3)
plot(tout,Y)
legend('x', 'y', 'z');
xlabel('t')
title(solve)
subplot(2,2,4)
plot3(Y(3,:), Y(2,:), Y(1,:))
xlabel('z')
ylabel('y')
zlabel('x')
title(solve)

err = max(max(abs(Y-Ym')));
fprintf('\nh = %1.1e \t error = %.10e\n\n', h, err);

