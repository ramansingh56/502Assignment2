function [a,e,i,omega,w,f] = rv2oe(r_vec,v_vec,mu)

    % Advanced Orbital Mechanics - convert position & velocity to orbital elements

    i_vec = [1 0 0];
    j_vec = [0 1 0];
    k_vec = [0 0 1];
    r = norm(r_vec);
    v = norm(v_vec);
    
    % vis-viva equation for semi-major axis (a)

    a = 1/((2/r)-((v^2)/mu));

    % eccentricity

    e_vec = (((v^2)/mu)-(1/r))*r_vec - ((1/mu)*(dot(r_vec,v_vec))*v_vec);
    e = norm(e_vec);

    % angular momentum

    h_vec = cross(r_vec,v_vec);
    h = norm(h_vec);

    % inclination
    
    i = acos(dot((h_vec/h),k_vec));

    % node vector

    n_vec = cross(k_vec,h_vec);
    n = norm(n_vec);

    % Longitude of the ascending node

    omega = acos((dot(n_vec,i_vec))/n);

    if dot(n_vec,j_vec) < 0
        omega = 2*pi - omega;
    end

    % argument of periapse

    w = acos((dot(n_vec,e_vec))/(n*e));

    if dot(e_vec,k_vec) < 0
        w = 2*pi - w;
    end

    % true anomaly

    f = acos((dot(r_vec,e_vec))/(r*e));

    if dot(r_vec,v_vec) < 0
        f = 2*pi - f;
    end

end