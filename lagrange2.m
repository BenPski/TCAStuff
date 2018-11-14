function out = lagrange2(t,x,dT)
    delta = x(1);
    delta_dot = x(2);
   
    l = 0.248;
    L = 0.055;
    theta0 = 2*pi*188;
    N = 2;
    D = 2e-3-0.55e-3;
    d = 0.55e-3;
    phi = l*sqrt(1-L^2/l^2)/(D/2);
    
    J = pi*d^4/32;
    g = 9.81;
    mu = 1e9;
    E = 600e6;
    rho = 4e-4;
    G = E/3;
    I = J/2;
    
    m = 20e-3;
    
    delta_gamma = phi*delta/l^2-theta0/l*(1-1/(1+rho*dT(t)));
    
    dU = l/2*(G*J*2*delta_gamma*phi/l^2 + E*I*phi^2/l^4*(-2*(L+delta)*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2))/(sqrt(l^2-(delta+L)^2))));
    
    delta_ddot = g-dU/m-mu*J*phi/(m*l^2)*delta_dot;
    
    
    out = [delta_dot;delta_ddot];
end