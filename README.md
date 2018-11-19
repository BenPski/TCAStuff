# TCA Stuff

Some of the code associated with the TCA behavior. Includes both Dynamic and Static models.


For the dynamics it is just the system of ODEs representing the TCA.

To use it the geometry parameters can be editted in the files (should make the parameters aguments) and then use one of the ode solvers available in matlab. Also, switching between extension and contraction is done by changing the sign of the strain due to temperature, this equates to changing the sign of theta0.

Since I haven't taken the time to generalize the system entirely (mostly due to directional gravitational effects) there are a couple functions to use. 

Hanging experiment:
```matlab
sol = ode45(@(t,x) tcaHanging(t,x,@(t) 0, @(t) 0, 20e-3), [0, T], [25;0;0])
res = deval(sol,linspace(0,T,1000));
plot(linspace(0,T,1000),res(2,:))
```
The system variables are y=(temperature, delta, delta_dot). The function looks like f(t,x,F(t),V(t),M) with t as the time, x as the state variables, F as the force applied to the TCA (axially), V as the voltage applied, and M as the hanging mass.

The tca dynamics for somthing like when it is embedded in the body (has no hanging mass) is:
```matlab
sol = ode45(@(t,x) fullTCADynamics(t,x,@(t) 0, @(t) 0.5*(sin(20*t/T)+1)), [0, T], [25;0;0])
```
The arguments are the same as before except there is no mass. Importantly this function doesn't properly handle a hanging weight, it would have to be modified to have F be a function like F(x_ddot,x_dot,t), which messes things up since it needs to solve the equation for x_ddot. So, this is mostly useful for later.
