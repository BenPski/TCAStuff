function out = tcaDynamicsAgain(t,x,F,T,V)
    %attempting to figure out the full dynamics including the twisting
    %angle phi
    %incorporates the thermal and dynamics model
    
    temp = x(1);
    delta = x(2);
    delta_dot = x(3);
    phi = x(4);
    phi_dot = x(5);
    
    
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
    m_delta = density*Vol/4
    alpha = asin(delta/l+L/l);
%     m_phi = density*L*(pi/4*(2*d^4+d^2*(D/2)^2)+2*pi*cos(alpha)*((D/2+d/2)^4-(D/2-d/2)^4));
    r = l*cos(alpha)/phi;
    m_phi = density*l^2*d^2*cos(asin(L/l))*(4*r^2+3*d^2)/(32*r);
    
    dT = 1/(density*Cp)*(V(t)^2/(R*l*Vol)-h*(temp-Tamb));
    
    delta_tau = phi*delta/l^2-theta0*(1-1/(1+rho*(temp-Tamb)))/l;
    delta_kappa = phi/l^2*(sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2));
    
    dTdD = phi/l^2;
    dTdP = delta/l^2;
    dKdD = -(delta+L)/sqrt(l^2-(delta+L)^2);
    dKdP = (sqrt(l^2-(delta+L)^2)-sqrt(l^2-L^2))/l^2;
    
    dUdD = l*(G*J*delta_tau*dTdD+E*I*delta_kappa*dKdD);
    dUdP = l*(G*J*delta_tau*dTdP+E*I*delta_kappa*dKdP);
    
    dDTdD = phi/l^2;
    dDTdP = delta/l^2;
        
%     if delta == 0
%         phi_dot = 0;
%     else
%         phi_dot = (T(t)-dUdP)/(2*l*mu*J*dDTdP-phi*delta_dot/l^2+theta0*rho*dT/(l*(1+rho*(temp-Tamb))^2))*l^2/delta;
%     end
    
    delta_tau_dot = (phi_dot*delta+phi*delta_dot)/l^2-theta0*rho*dT/(l*(1+rho*(temp-Tamb))^2);
    
    dEdD = 2*l*mu*J*delta_tau_dot*dDTdD;
    dEdP = 2*l*mu*J*delta_tau_dot*dDTdP;
    
   
    delta_ddot = (F(t)-dUdD+m_delta*g-dEdD)/m_delta;
    phi_ddot = (T(t)-dUdP-dEdP)/m_phi;
    
    
    
    out = [dT;delta_dot;delta_ddot;phi_dot;phi_ddot];
%     out = [dT;delta_dot;delta_ddot;phi_dot];
end