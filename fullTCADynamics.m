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
    dT = temp-Tamb;
    l = 0.248;
    d = 0.55e-3;
    vol = l*pi*d^2/4;
    
    temp_dot = (V(t)^2/(R*l*vol)-h*(temp-Tamb))/(density*Cp);
    
    L = 0.055;
    theta0 = 2*pi*188;
    N = 2;
    D = 2e-3-0.55e-3;
    phi = l*sqrt(1-L^2/l^2)/(D/2);
    
    J = pi*d^4/32;
    mu = 1e9;
    E = 600e6;
    rho = 4e-4;
    G = E/3;
    I = J/2;
    
    
    delta_gamma = phi*delta/l^2-theta0/l*(1-1/(1+rho*dT));
    
    dU = l/2*(G*J*2*delta_gamma*phi/l^2 + E*I*phi^2/l^4*(-2*(L+delta)*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2))/(sqrt(l^2-(delta+L)^2))));
    
    delta_ddot = (4/(density*vol))*(F(t)-dU-mu*J*phi/l^2*delta_dot);
    
    out = [temp_dot;delta_dot;delta_ddot];
    
end