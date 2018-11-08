function out = lagrangeDynamics(t,d,dT,mu)
    %the lagrange dynamics formulation
    
    delta = d(1);
    delta_dot = d(2);
   
    l = 0.248;
    L = 0.055;
    phi0 = 2*pi*188;
    N = 2;
    D = 2e-3-0.55e-3;
    d = 0.55e-3;
    
    E = 600e6;
    nu = 0.45;
    rho = (4.2e-04);
    %mu = 1e9;
    
    d0 = d/(2+sin(pi/N));
    
    alpha0 = asin(L/l);
    
    d = d0*(1+rho*dT);
    phi = d0*phi0/d;
    gamma = (phi0-phi)/l;
%     alpha = asin(delta/l+sin(alpha0));
%     gamma = 2*cos(alpha0)/D*(1/tan(alpha)-1/tan(alpha0));
    G = E/(2*(1+nu));
        
    A = N*pi/4*d^2;
    J = N*pi/4*(sin(pi/N)^2+1/8)*d^4;
    I = J/2;
    n = l*cos(alpha0)/(pi*D);
    
    A1 = 8*n*(l/n)^3/(pi^3*d^4*G);
    A2 = 8*n*(l/n)/(pi*d^2*2*G);
    A3 = 8*n*(l/n)^3*2/(pi^3*d^4*E);
    A4 = 8*n*(l/n)/(pi*d^2*2*E);
    B1 = 8*n*(l/n)^2/(pi^2*d^4*G);
    
    d_f12 = -(2*B1*(sin(alpha0) + delta/l))/l;
    d_f11 = -(2*(delta + l*sin(alpha0))*(2*A3*delta^2 - 2*A1*delta^2 + A1*l^2 + A2*l^2 - A4*l^2 + A1*l^2*cos(2*alpha0) - A3*l^2*cos(2*alpha0) - 4*A1*delta*l*sin(alpha0) + 4*A3*delta*l*sin(alpha0)))/l^4;
    
    f11 = A1*((sin(alpha0) + delta/l)^2 - 1)^2 - A2*((sin(alpha0) + delta/l)^2 - 1) + A4*(sin(alpha0) + delta/l)^2 - A3*((sin(alpha0) + delta/l)^2 - 1)*(sin(alpha0) + delta/l)^2;
    f12 = -B1*((sin(alpha0) + delta/l)^2 - 1);
    f22 = 32*l/(pi*d^4*G);
    
    alpha = asin(delta/l+sin(alpha0));
    Q_nc = -mu*J*cos(alpha)/(l^2*D/2)*delta_dot;
    g = -9.81;
    m = 10e-3;
    delta;

    %gamma = cos(alpha)/(D/2)*(delta/l+sin(alpha0));
    tau = G*J*gamma;
    
    U = 1/2*f11*(m*g)^2-2*f12*(m*g)*tau+1/2*f22*tau^2
    
    %delta_ddot = g + 1/2*d_f11*m*g^2 - 2*d_f12*g*tau + Q_nc/m;
    delta_ddot = fsolve(@(delta_ddot) -g + 1/2*d_f11*(m*(g-delta_ddot))^2/m - 2*d_f12*(g-delta_ddot)*tau + Q_nc/m - delta_ddot,0);
    
    
    out = [delta_dot;delta_ddot];
end