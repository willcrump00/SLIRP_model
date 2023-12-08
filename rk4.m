function [y] = rk4(odefun,t,dt,y0,e,mu_L,p)
%calculating partial slopes
k1 = odefun(t,y0,e,mu_L,p);
k2 = odefun(t  ,y0 +((.5*k1)*dt),   e,mu_L,p);
k3 = odefun(t  ,y0 +((.5*k2)*dt),   e,mu_L,p);
k4 = odefun(t  ,y0 +(k3*dt),   e,mu_L,p);
%combining to calulate final integrated value of y
y = y0 + ((1/6) * (k1+2*k2+2*k3+k4)*dt);
end