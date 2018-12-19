% Initialization

clc
close all
clear

format long g

% Data
mu_d = 398600.;       % Earth gravitational parameter [km^3/s^2]
Re_d = 6378.;         % Earth radius [km]

% Conversion units
rconv = Re_d;
vconv = sqrt(mu_d / Re_d);
tconv = rconv / vconv;
mu = 1.;

% Initial conditions
x0 = (Re_d + 700.) / rconv;
y0 = 0.;
z0 = 0.;
vx0 = 0.;
vy0 = sqrt(mu / x0);
vz0 = 0.;
X0 = [x0, y0, z0, vx0, vy0, vz0]';

% Time
t0 = 0.;
t1 = 8000. / tconv;
tspan = [t0, t1];

% Ode fun
ode_fun = @(t, y) ode_kepler(t, y, mu);

% Options
AT = 1e-3;
RT = 1e-3;
options = odeset('AbsTol', AT, 'RelTol', RT);

% Numerical integration
[t, MAT] = ode45(ode_fun, tspan, X0, options);

% Print results
filename = 'mat_results.dat';
MAT = [t MAT];

print_results(MAT, filename);
