% Advanced Orbital Mechanics Assignment 2 Problem 1 & 2
% Raman Singh
% Molniya orbit for Earth & Mars
% All computations are done in metric units

close all; clear; clc;

% choose 1 for Earth and 2 for Mars
planet = 2;

switch planet
    case 1
        % for Earth
        mu = 3.986e5;
        R = 6370;
        J2 = 0.00108;
        day = 3600*24;
        P = 8*3600;
        rp = 600;

    case 2
        % for Mars
        mu = 4.282e4;
        R = 3390;
        J2 = 0.00196;
        P = 24*3600 + 39*60 + 35;
        day = P;
        rp = 400;
end

% semi-major axis
a = ((P*sqrt(mu))/(2*pi))^(2/3);

% mean motion
n = sqrt(mu/a^3);

% eccentricity (varying from 0 to its value corresponding to the periapse
% limit specified in the problem)
e = linspace(0,1 - ((rp+R)/a),5000);

% solution of rate of change of avg argument of periapse
inc = [acos(1/sqrt(5)) acos(-1/sqrt(5))];

for i = 1:length(inc)
    for j = 1:length(e)
        % nodal precession drift rate
        omega_dot(i,j) = (-1.5*J2*n*((R/a)^2)*(cos(inc(i))))/((1 - e(j)^2)^2);
        if i == 1
            % velocity at apoapse
            vel_apo(j) = sqrt((mu/a)*((1-e(j))/(1+e(j))));
        end
    end
end

% plotting
figure
plot(e,omega_dot(1,:)*180/pi,LineWidth=1.5)
xlabel("Eccentricity (e)")
ylabel('Nodal Precession Drift Rate $\dot{\overline{\Omega}}$ in deg/s','Interpreter','latex')
title("Variation of Nodal Precession Drift Rate","with eccentricity")

figure
plot(e,omega_dot(1,:)*day*180/pi,LineWidth=1.5)
xlabel("Eccentricity (e)")
ylabel('Nodal Precession Drift Rate $\dot{\overline{\Omega}}$ in deg/day','Interpreter','latex')
title("Variation of Nodal Precession Drift Rate with eccentricity")

figure
plot(a*(1-e),omega_dot(1,:)*180/pi,LineWidth=1.5)
xlabel("Radius of periapsis (r_p) in km")
ylabel('Nodal Precession Drift Rate $\dot{\overline{\Omega}}$ in deg/s','Interpreter','latex')
title("Variation of Nodal Precession Drift Rate","with radius of periapsis")

figure
plot(e,vel_apo,LineWidth=1.5)
xlabel("Eccentricity (e)")
ylabel("Velocity at apoapsis (v_a) in km/s")
title("Variation of eccentricity with velocity at apoapsis")