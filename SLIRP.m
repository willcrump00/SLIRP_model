function [s,l,I,r,p] = SLIRP(params, tspan , T )
global ui;
global er;
global ap;
global binst;
global instantdpb;
global ul;
global instantdp;
instantdp = 0;
betamax = params(1);
ui = params(2);
er = params(3);
ap = params(4);
ul = params(5);
%preallocating for speed (GOTTA GO FAST)
n = length(tspan);
p = zeros(1,n);
s = zeros(1,n);
l = zeros(1,n);
I = zeros(1,n);
r = zeros(1,n);
pb = zeros(1,n);
%CALCULATE INITIAL VALUES
p(1) = 1.33*30*(-0.35968 + (.10789*15)-.00214*15*15)*30;
s(1) = p(1)/ap;
l(1) = 0.01*s(1);
I(1) = 0;
r(1) = ui*I(1);
b(1) = 1;

%find Te
Te = -.35968 + (.10789.*T) - (.00213*(T.^2));

%berries
pb(1) = (.1724*b(1)-0.0000212*b(1)^2)*Te(1);

%find tb
Tb = zeros(1,length(T));
for i = 1:length(tspan)
    if T(i)>0 && T(i)<35
            Tb(i) = 0.000241*T(i)^2.06737*(35-T(i))^.72859;
    else
        Tb(i) = 0;
    end
end
%find beta & ul-1
beta = betamax .*Tb;
binst = beta(1);

% dpldt is unchanged by any variables so we can calc outright
dpl = (1.33.*tspan).*Te;

%equations
%just changed pb
dpb = @(pb,Te) (0.1724*pb - (.0000212*(pb^2))) *Te;
dp = @(p, instantdpb, dpl) instantdpb+dpl;

%SLIRP EQUATIONS
ds = @(s,I) (-1*binst)*(s*I) + (instantdp * (1/ap));
dl = @(l,I,s) (s*I) - ((1/ul)*l) + er;
di = @(I,l) ((1/ul)*l) - ((1/ui)*I);
dr = @(r, I) (1/ui)*I;

%slope calculation
for i = 1:n-1
    pb(i+1) = rk4(dpb,(tspan(i+1)-tspan(i)),[pb(i),Te(i)]);
    %ghetto problems require ghetto solutions
    %introducing! the global variable catastrophe
    instantdpb = dpb(pb(i),Te(i));
    binst = beta(i);
    p(i+1) = rk4(dp,(tspan(i+1)-tspan(i)),[p(i),instantdpb, dpl(i)]);
    instantdp = dp(p(i),instantdpb,dpl(i));
    s(i+1) = rk4(ds, (tspan(i+1)-tspan(i)), [s(i),I(i)]);
    l(i+1) = rk4(dl, (tspan(i+1)-tspan(i)), [l(i),I(i),s(i)]);
    I(i+1) = rk4(di, (tspan(i+1)-tspan(i)), [I(i), l(i)]);
    r(i+1) = rk4(dr,(tspan(i+1)-tspan(i)), [r(i), I(i)]);
end
end

