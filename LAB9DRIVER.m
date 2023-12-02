%
%DP/Dt = k(T,t)
%T temp t days since event
%Te = -.35968 + .10789T - .00213T^2
% T air temp - - Te 0-1 growth rates
%Dpb/dt (172.4Pb - 21.2Pb^2)*Te
%pb = berry surface area
%DPl/dt - (1.33tday)*Te
%dpl is leaf surface area
%dp/dt = dpb/dt + dpl/dt
%Tb = 0 if T<=0
%Tb = .000241T^2.06737(35-T)^.72859 0<T<35
%Tb = 0 if T>=35
%beta = betamax*Tb;
%ul = sum Tb from tinitial to tend
clc
clear
close all

load EnvironmentalForcing.mat

bmax = 1;
ulmin = 6;
ui = 10;
er = 0.001;
ap = 5000;
p = 1.33*30*(-0.35968 + (.10789*15)-.00214*15*15)*30;
si = p/ap;
li = 0.01*si;
ii = 0;
ri = ui*ii;
bi = 0;


%dsdt = -beta*s*i + dpdt*(1/ap);
%dldt = s*i - (ul*l+er);
%didt = (ul*l) - (ui*i);
%dpdt = dpbdt + dpldt;
%dpbdt = (172.4Pb - 21.2Pb^2)*Te;

% input pararms are : betamax, ui, er, ap, ulmin
params = [bmax,ui,er,ap,ulmin];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');
%copy tb calcs
Tb = zeros(1,length(T));
for i = 1:length(tspan)
    if T(i)>0 && T(i)<35
            Tb(i) = 0.000241*T(i)^2.06737*(35-T(i))^.72859;
    else
        Tb(i) = 0;
    end
end
figure
plot(T,Tb,'*');
xlabel('air temp [C]');
ylabel('Tb');
title('effect of air temp on Tb');

Te = -.35968 + (.10789.*T) - (.00213*(T.^2));
figure
plot(tspan,Te,'-');
xlabel('time [days]');
ylabel('Tb');
title('Effect of temperature on forcing air temp');

%running different b max
params = [10,ui,er,ap,ulmin];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model - bmax = 10');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');

params = [.1,ui,er,ap,ulmin];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model - bmax = .1');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');

%3 different different ul min
params = [bmax,ui,er,ap,2];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model ulmin = 2');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');

params = [bmax,ui,er,ap,15];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model - ulmin = 15');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');

params = [bmax,ui,er,ap,0.1];
[s,l,i,r,p] = SLIRP(params, tspan, T);
p= p./ap;
figure
plot(tspan,s,tspan,l,tspan,i,tspan,r,tspan,p);
title('plant disease model - ulmin = 0.1');
ylabel('fraction of population');
xlabel('time [days]');
legend('green: pop , blue: S , red: L , orange: I , purple: R','Location', 'northwest');


