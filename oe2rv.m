function [r,r_vec,v,v_vec] = oe2rv(a,e,i,omega,w,M,mu)

    % Advanced Orbital Mechanics - convert orbital elements to position and velocity

    % basis vectors
    I_vec = [1 0 0];
    J_vec = [0 1 0];
    K_vec = [0 0 1];

    % convert the angle from degree to radians
    i = i * (pi/180);
    omega = omega * (pi/180);
    w = w * (pi/180);
    M = M * (pi/180);

    % finding true anomaly using mean anomaly via eccentric anomaly
    E = M;
    tol = 1;
    while tol > 1e-15
        Enew = E - (E - e*sin(E) - M)/(1 - (e*cos(E)));
        tol = Enew - E;
        E = Enew;
    end

    E_out = E;
    f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

    % angle of rotation
    th = w + f;

    % magnitude of the radius vector
    r = (a*(1-(e^2)))/(1+(e*cos(f)));

    % position vector
    r_vec = r*((cos(th)*cos(omega) - cos(i)*sin(omega)*sin(th))*I_vec + (cos(th)*sin(omega) + cos(i)*cos(omega)*sin(th))*J_vec + (sin(i)*sin(th))*K_vec);

    % angular momentum
    h = sqrt(mu*r*(1 + e*cos(f)));

    % velocity vector
    v_vec = - (mu/h)*((cos(omega)*(sin(th)+e*sin(w))) + (sin(omega)*(cos(th)+e*cos(w))*cos(i)))*I_vec ...
        - (mu/h)*((sin(omega)*(sin(th)+e*sin(w))) - (cos(omega)*(cos(th)+e*cos(w))*cos(i)))*J_vec ...
        + (mu/h)*(cos(th)+e*cos(w))*sin(i)*K_vec;

    % magnitude of the velocity vector
    v = norm(v_vec);
    
end