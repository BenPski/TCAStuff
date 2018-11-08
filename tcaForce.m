function out = tcaForce(delta,temp)
    %alright this time I think the model should be correct
    %the last fix essentially assumed the diameter of the helix was
    %constant which is wrong and I think the initial one may have some
    %mistake in the terms
    %however no analytic solution for (F,temp) -> delta, have to
    %iteratively invert (delta,temp) -> F
    
    l = 0.248;
    L = 0.055;
    phi0 = 2*pi*188;
    N = 2;
    D = 2e-3-0.55e-3;
    d = 0.55e-3;
    
    %E = (600-0.4655*(temp+25)^1.381)*10^7;
    E = (600-4*temp)*10^7;
    nu = 0.45;
    rho = (4.2e-04+5e-6*(temp+25));
    
    d0 = d/(2+sin(pi/N));
    
    alpha0 = asin(L/l);
    
    d = d0*(1+rho*temp);
    phi = d0*phi0/d;
    gamma = (phi0-phi)/l;
%     alpha = asin(delta/l+sin(alpha0));
%     gamma = 2*cos(alpha0)/D*(1/tan(alpha)-1/tan(alpha0));
    G = E/(2*(1+nu));
        
    A = N*pi/4*d^2;
    J = N*pi/4*(sin(pi/N)^2+1/8)*d^4;
    I = J/2;
    
    tau = G*J*gamma;
    
    n = l*cos(alpha0)/(pi*D);
    alpha = asin(delta/l+sin(alpha0));
    f1 = 1/2*l*(l^2/(G*J*4*pi^2*n^2)*cos(alpha)^4 + l^2/(E*I*4*pi^2*n^2)*cos(alpha)^2*sin(alpha)^2 + cos(alpha)^2/(G*A) + sin(alpha)^2/(E*A));
    f2 = l^2/(G*J*2*pi*n)*cos(alpha)^2;
    F = (delta+tau*f2)/f1;
    out = F;
    