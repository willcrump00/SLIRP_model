function [y] = rk4(fun,dx,y0,vars)
%figure out input length for different functin
check = length(y0);
%only iterates once per calling, so no loop needed
if check == 1
k1 = fun(y0(1));
k2 = fun(y0(1)+(.5*k1)*dx);
k3 = fun(y0(1)+(.5*k2)*dx);
k4 = fun(y0(1)+k3*dx);
else if check ==2
k1 = fun(y0(1),y0(2));
k2 = fun(y0(1)+(.5*k1)*dx , y0(2)+(.5*k1)*dx);
k3 = fun(y0(1)+(.5*k2)*dx , y0(2)+(.5*k2)*dx);
k4 = fun(y0(1)+k3*dx , y0(2)+k3*dx);
else
k1 = fun(y0(1),y0(2),y0(3));
k2 = fun(y0(1)+(.5*k1)*dx , y0(2)+(.5*k1)*dx , y0(3)+(.5*k1)*dx);
k3 = fun(y0(1)+(.5*k2)*dx , y0(2)+(.5*k2)*dx , y0(3)+(.5*k2)*dx);
k4 = fun(y0(1)+k3*dx , y0(2)+k3*dx , y0(3)+k3*dx);
end

y = y0(1) + (1/6) * (k1+2*k2+2*k3+k4)*dx;
end