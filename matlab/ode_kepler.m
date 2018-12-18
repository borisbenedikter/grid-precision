function dX = ode_kepler(t, X, mu)
%ode_kepler - ODEs for a Keplerian orbit.
%
% Syntax: dX = ode_kepler(t, X, mu)
%

% Input
x = X(1);
y = X(2);
z = X(3);
vx = X(4);
vy = X(5);
vz = X(6);

r = sqrt(x^2 + y^2 + z^2);

% Output
dX = zeros(length(X), 1);
dX(1) = vx;
dX(2) = vy;
dX(3) = vz;
dX(4) = -mu / r^3 * x;
dX(5) = -mu / r^3 * y;
dX(6) = -mu / r^3 * z;

end