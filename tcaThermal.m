function out = tcaThermal(t,T,V)
    %the thermal ode for a tca
    
    rho = 1300; %density
    Cp = 1267.7; %specific heat
    R = 3.5; %resistance per length
    h = 20; %convection with air
    Tamb = 25;
    l = 0.248;
    d = 0.55e-3;
    vol = l*pi*d^2/4;
    
    out = (V(t)^2/(R*l*vol)-h*(T-Tamb))/(rho*Cp);
end