function out = fullerTCADynamics(t,x,F,T,V)
    %compute the full dynamics of the TCA
    %represents the ode for the displacement oand the temperature
    %V(t) is input voltage, F(t) is the applied force (positive is in
    %direction of positive displacement)
    
    %also includes phi dynamics (use a change to phi to make it easier to
    %initialize)
    
    temp = x(1);
    delta = x(2);
    delta_dot = x(3);
    dphi = x(4);
    dphi_dot = x(5);
    
    density = 1300; %density
    Cp = 1267.7; %specific heat
    R = 3.5; %resistance per length
    h = 20; %convection with air
    Tamb = 25;
    dT = temp-Tamb;
    l = 0.248;
    d = 0.55e-3;
    vol = l*pi*d^2/4;
    
    temp_dot = (V(t)^2/(R*l*vol)-h*(temp-Tamb))/(density*Cp);
    
    L = 0.055;
    theta0 = 2*pi*188;
    N = 2;
    D = 2e-3-0.55e-3;
    phi0 = l*sqrt(1-L^2/l^2)/(D/2); %the initial phi
    phi = phi0+dphi;
    
    J = pi*d^4/32;
    mu = 1e6;
    %mu = 0.1*l^2/(phi*J)
    E = 600e6;
    rho = 4e-4;
    G = E/3;
    I = J/2;
    
    
    delta_gamma = phi*delta/l^2-theta0/l*(1-1/(1+rho*dT));
    delta_kappa = phi/l^2*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2));
    
    %dUdD = l/2*(G*J*2*delta_gamma*phi/l^2 + E*I*phi^2/l^4*(-2*(L+delta)*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2))/(sqrt(l^2-(delta+L)^2))));
    dUdD = l/2*(G*J*2*delta_gamma*phi/l^2 + 2*E*I*delta_kappa*phi/l^2*(-(delta+L)/sqrt(l^2-(delta+L)^2)));
    
    delta_ddot = (4/(density*vol))*(F(t)-dUdD-mu*J/l^3*(2*phi^2*delta_dot+2*phi*delta*dphi_dot));
    %delta_ddot = (4/(density*vol))*(F(t)-dUdD);
    
    dUdP = l/2*(2*G*J*delta_gamma*delta/l^2 + 2*E*I*delta_kappa*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2))/l^2);
    
    dphi_ddot = (T(t)-dUdP-mu*J/l^3*(2*delta^2*dphi_dot+2*phi*delta*delta_dot))/(density*L*pi/32*(D^4-(D-d)^4));
    
    out = [temp_dot;delta_dot;delta_ddot;dphi_dot;dphi_ddot];
    
end