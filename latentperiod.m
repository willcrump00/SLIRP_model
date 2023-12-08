% latent period length as a function of Temperature (Calonnec et al. 2008)
function mu_L = latentperiod(istart,dt,Nsteps,mu_L_target,mu_L,T)


%latent period length as a function of Temperature (Calonnec et al. 2008)
%calculate temp effect array.  This could be calculated once external to
%time loop for efficiency and passed to this function
PT   = zeros(size(T));
for i=1:Nsteps
    PT(i)=Sall_temp_effect(T(i));
end
%now calculate the time to latent from istart (when a new vine was infected)
flag=true;
for i=istart:Nsteps
    if(flag)
        PTint = trapz(PT(istart:i))*dt;
        if(PTint >= mu_L_target)
            mu_L(istart:i)=dt*(i-istart);
            ispan = i;
            istart= istart+1;
            flag  = false;
        end
    else
        PTint = trapz(PT(istart:i))*dt;
        if(PTint >= mu_L_target)
            mu_L(ispan:i)=dt*(i-istart);
            istart=istart+1;
            ispan=i;
        end
    end
end
mu_L=1./mu_L;

%set any infinite mu_L's back to zero
infInd=isinf(mu_L);
mu_L(infInd)=0;

end