%% Advanced Orbital Mechanics Assignment 2 Problem 3
% Raman Singh
% Special Perturbation Method -- Cowell's Method
% Inspired from Dr. Woollands lecture on J2 perturbations & Curtis
% All computations are done in metric units

close all; clear; clc;

%% Initialization
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode

% Orbital elements
a = 26600;
e = 0.74;
i = 1.10654*180/pi;
omega = 90;
w = 5;
M = 10;

% Gravitational parameter
mu = 398600;

[~,r_0,~,v_0] = oe2rv(a,e,i,omega,w,M,mu);

% format long
% disp(r_0);
% disp(v_0);

x0 = [r_0 v_0];

P = 2*pi*sqrt((a^3)/mu);

t0 = 0;
tf_days = 100;
tf = tf_days*3600*24;
tspan = linspace(t0,tf,50000);
day = 3600*24;

% numerical integration
[t,v] = ode45(@(t,v) two_bodyJ2(t,v), tspan, x0, opts_ode);

for j = 1:length(v)
    [a,e,i,omega,w,f] = rv2oe([v(j,1:3)],[v(j,4:6)],mu);
    semim(j) = a;
    ecc(j) = e;
    inc(j) = i;
    lan(j) = omega;
    ap(j) = w;
    ta(j) = f;
end

% plotting
figure
plot3(v(:,1),v(:,2),v(:,3),LineWidth=0.5,Color='r');
hold on
grid on
scatter3(0,0,0,1000,"blue",'filled');
xlabel('x in km')
ylabel('y in km')
zlabel('z in km')
title('Orbit Plot')

% variation of orbital elements with time
figure
subplot(3,2,1)
plot(t/day,semim,LineWidth=1);
xlabel('t in days')
ylabel('Semimajor axis (a) in km')
title('Variation of semi-major axis with time')

subplot(3,2,2)
plot(t/day,ecc,LineWidth=1);
xlabel('t in days')
ylabel('Eccentricity (e)')
title('Variation of eccentricity with time')

subplot(3,2,3)
plot(t/day,inc*(180/pi),LineWidth=1);
xlabel('t in days')
ylabel('Inclination (i) in ^o')
title('Variation of inclination with time')

subplot(3,2,4)
plot(t/day,lan*(180/pi),LineWidth=1);
xlabel('t in days')
ylabel('Right ascension of the ascending node (\Omega) in ^o')
title('Variation of right ascension of the ascending node with time')

subplot(3,2,5)
plot(t/day,ap*(180/pi),LineWidth=1);
xlabel('t in days')
ylabel('Argument of periapsis (\omega) in ^o')
title('Variation of argument of periapsis with time')

subplot(3,2,6)
plot(t/day,ta*(180/pi),LineWidth=1);
xlabel('t in days')
ylabel('True anomaly (f) in ^o')
title('Variation of true anomaly with time')

%% Two body problem differential equations with J2 Perturbations

function dvdt = two_bodyJ2(t,v)

    mu = 398600;
    R = 6370;
    J2 = 0.00108;

    r = sqrt((v(1))^2+(v(2))^2+(v(3))^2);
    dvdt = zeros(6,1);
    
    % J2 perturbation
    k = (1.5*J2*mu*(R^2))/(r^4);

    p_x = k*(v(1)/r)*((5*(v(3)/r)^2)-1);
    p_y = k*(v(2)/r)*((5*(v(3)/r)^2)-1);
    p_z = k*(v(3)/r)*((5*(v(3)/r)^2)-3);

    % differentiation of the x distance wrt time
    dvdt(1) = v(4);

    % differentiation of the y distance wrt time
    dvdt(2) = v(5);

    % differentiation of the z distance wrt time
    dvdt(3) = v(6);

    % differentiation of the x velocity wrt time
    dvdt(4) = - (mu*v(1))/(r^3) + p_x;

    % differentiation of the y velocity wrt time
    dvdt(5) = - (mu*v(2))/(r^3) + p_y;

    % differentiation of the z velocity wrt time
    dvdt(6) = - (mu*v(3))/(r^3) + p_z;
    
end