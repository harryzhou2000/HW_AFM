

FR = @(t,f)  2 * (f - f^2);
options = odeset('InitialStep',1e-5,'RelTol',1e-8,'AbsTol',1e-8);
sol = ode45(FR, [0,10], 0.000000001, options);
plot(sol.x, sol.y,'-o');

%%
syms t
syms f(t)
S = dsolve(diff(f,t) == 2 * (f-f^2), f(0) == 1/3)

% -1/(exp(C4 - 2*t) - 1)
