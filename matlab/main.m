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

% Initial conditions [non-dim]
a_orb = (Re_d + 700.) / rconv;
e_orb = 0.1;
i_orb = deg2rad(45);
gom = deg2rad(30);
pom = deg2rad(30);
nu = deg2rad(15);
COE = [a_orb, e_orb, i_orb, gom, pom, nu];
[rECI0, vECI0] = coe2rvECI(COE, mu);

x0 = rECI0(1);
y0 = rECI0(2);
z0 = rECI0(3);
vx0 = vECI0(1);
vy0 = vECI0(2);
vz0 = vECI0(3);

% x0 = (Re_d + 700.) / rconv;
% y0 = 0.;
% z0 = 0.;
% vx0 = 0.;
% vy0 = sqrt(mu / x0);
% vz0 = 0.;

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
