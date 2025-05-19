function dy = fullODE(t, y)
    
    J = diag([2.5, 2.5, 3.3]);
    invJ = inv(J);
    r = [-0.001; 0; -0.001];
    m = 55;
    g=-9.81;
   

   
    Omega = y(1:3);
    R = reshape(y(4:12), 3, 3);  % R is 3x3

    if abs(R(3,1)) < 1
        pitch = asin(-R(3,1));              
        yaw   = atan2(R(2,1), R(1,1));      
        roll  = atan2(R(3,2), R(3,3));       
    else
        pitch = pi/2 * sign(-R(3,1));
        yaw   = atan2(-R(1,2), R(2,2));
        roll  = 0;
    end
    E = [yaw; pitch; roll];

    g_vec = [-g * sin(E(2));
              g * sin(E(1)) * cos(E(2));
              g * cos(E(1)) * cos(E(2))];


 
    Omega_cross = [
        0        -Omega(3)  Omega(2);
        Omega(3)  0         -Omega(1);
        -Omega(2) Omega(1)   0
    ];
    coriolis = Omega_cross * (J * Omega);
    gravity_moment = cross(r, m * g_vec);
    dOmega = invJ * ( - gravity_moment - coriolis);
 

    Omega_hat = Omega_cross;
    dR = R * Omega_hat;
    
    dR_vec = reshape(dR, 9, 1);
   
    dy = [dOmega; dR_vec];
end