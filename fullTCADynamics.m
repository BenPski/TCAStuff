function out = fullTCADynamics(t,x,F,V)
    %compute the full dynamics of the TCA
    %represents the ode for the displacement oand the temperature
    %V(t) is input voltage, F(t) is the applied force (positive is in
    %direction of positive displacement)
    
    temp = x(1);
    delta = x(2);
    delta_dot = x(3);

    
    
    density = 1300; %density
    Cp = 1267.7; %specific heat
    R = 3.5; %resistance per length
    h = 20; %convection with air
    Tamb = 25;
    l = 0.248;
    d = 0.55e-3;
    Vol = l*pi*d^2/4;
    
    
    L = 0.055;
    theta0 = 2*pi*188;
    D = 2e-3-0.55e-3;

    
    J = pi*d^4/32;
    mu = 1e3;
    %mu = 0.1*l^2/(phi*J)
    E = 600e6;
    rho = 4e-4;
    G = E/3;
    I = J/2;
    
    g = 9.81;
    m_delta = density*Vol/4;
    phi = l*sqrt(1-L^2/l^2)/(D/2);
    
    
    dT = 1/(density*Cp)*(V(t)^2/(R*l*Vol)-h*(temp-Tamb));
    
    delta_tau = phi*delta/l^2-theta0*(1-1/(1+rho*(temp-Tamb)))/l;
    delta_kappa = phi/l^2*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2));
    
    dTdD = phi/l^2;
    dKdD = -(delta+L)/sqrt(l^2-(delta+L)^2);
    
    dUdD = l*(G*J*delta_tau*dTdD+E*I*delta_kappa*dKdD);
    
    dDTdD = phi/l^2;
    
    delta_tau_dot = (phi*delta_dot)/l^2-theta0*rho*dT/(l*(1+rho*(temp-Tamb))^2);
    
    dEdD = 2*l*mu*J*delta_tau_dot*dDTdD;
    
   
    delta_ddot = (F(t)-dUdD+m_delta*g-dEdD)/m_delta;
    
    out = [dT;delta_dot;delta_ddot];
    
end