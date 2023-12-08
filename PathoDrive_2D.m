% driver script for our 2 dimensional pathogen simulation for ME EN 2450
% Fall 2023.
%
close all;clear all;clc
%
% load Envinronmental Forcing data (U,V,T,tspan)
load EnvironmentalForcing.mat

%declare global varialbes so the function have access without passing
global NpX NpY Nsteps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set simulation constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Global domain constants (apply everywhere) %%%%%%%%%%%%%%%%%
NpX = 50;     %number of plants in the X-direction
NpY = 50;     %number of plants in the Y-direction
Nsteps = length(T); %number of time steps for integration

%%%%%%%%%%%%%%%%%%%%%%%%% Global plant parameters %%%%%%%%%%%%%%%%%%%%%%%%%
LpX = 1;      %X-direction physical length of a plant in meters
LpY = 1;      %Y-direction physical length of a plant in meters
%normalization factor for a plant ('final' susceptible plant surface area in cm^2)
A   = 5000;     
%initial average size of a individual plant in the population 
% (=model after 30 days with constant temp of 15C + 1 for initial berry size)
P_i_ave  = 1.33*30*(-0.35968+0.10789*15-0.00214*15*15)*30+1; 
P_i_std = 0.2*P_i_ave; %variance in the initial growth (fraction of the average)

%%%%%%%%%%%%%%%%%%%%%%% Global pathogen parameters %%%%%%%%%%%%%%%%%%%%%%%%
beta_max = 2;   %max rate infection spread under ideal conditions (1/day)
mu_L_min = 6;   %min length of latent period (min number of days latent)
mu_I  = 10;     %rate infection clears (number of days infectious)
eta   = 1;      %release fraction scale factor
kappa = 0.75;   %release fraction scale factor   
xi    = -2.0;   %release fraction offset
Gamma = 1e-2;   %spore production from Calonnec et al 2009 approx scaled as surface area coverage
alpha = 0.314;  %spore production 2nd factor

%%%%%%%%%%%%%%%%%% Initialize individual Plants (vines) %%%%%%%%%%%%%%%%%%%
% Here we will use a structure (vine) to store all the different variables
% to keep the association between variables and locations
vn = {'X','Y','IsInfect','LatentSwitch','P','B','S','L','I','R','E','F','mu_L'}; %list of variable names to store in the structure
% initialize the storage structure
vine(1:NpX*NpY) = cell2struct({0,0,false,false,zeros(Nsteps,1),zeros(Nsteps,1),...
    zeros(Nsteps,1),zeros(Nsteps,1),zeros(Nsteps,1),zeros(Nsteps,1),...
    zeros(Nsteps,1),zeros(Nsteps,1),zeros(Nsteps,1)},vn,2);
% set the position and initial variables.  Here we are using a vector for
% later convience when looping over only subsets of vines with active
% infections
for i=1:NpX
    for j=1:NpY
        cnt = i+(j-1)*NpX; %counter to vectorize vine 
        vine(cnt).X = i*LpX-0.5; %x-position of the center of a vine in meters
        vine(cnt).Y = j*LpY-0.5; %y-position of the center of a vine in meters
        vine(cnt).P(1) = P_i_ave+P_i_std*randn; %initial population size (random around the average)
        vine(cnt).B(1) = 1;        %initial size of the berry population (assumed small (1cm^2))
        vine(cnt).S(1) = vine(cnt).P(1)/A;    %initial size of the susceptible population (normalized)
        vine(cnt).L(1) = 0;        %initial fraction of the population that is latent
        vine(cnt).I(1) = 0;        %initial fraction of the population that is infectious 
        vine(cnt).R(1) = vine(cnt).I(1)*mu_I; %initial fraction of the population that is recovered
        vine(cnt).E(1) = 0;        %initial fraction of population from external sources
        vine(cnt).F(1) = 0;        %initial amount of 'spores' for spreading
    end
end

% Randomly select an initial plant for the infection to start at and give
% it the equivalent of 1 small (0.5cm diameter) spot on a leaf that is infected (latent)
RandV = randi(NpX*NpY); %a random integer value in the range of the number of vines
vine(RandV).L(1) = 0.25^2/4*pi/A;

% call the pathogen function
tic
[vine]=PathogenGrowth_2D(vine,beta_max,mu_L_min,mu_I,A,eta,kappa,xi,Gamma,...
    alpha,T,U,V,tspan);
toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Calculate stats for entire domain at each time step
S_ave = mean([vine.S],2);
%INSERT YOUR CODE HERE to fill in the rest...

%%%%%%%%%%%%%%%%%%%%% Plot the average of the field %%%%%%%%%%%%%%%%%%%%%%%
FSize = 14; %fontsize for plots
figure;plot(tspan,S_ave,'-k','LineWidth',2);
legend({'Susceptible'},'Location','NorthWest');
xlabel('time (days)','Fontsize',FSize);
ylabel('Population (fraction of initial)','Fontsize',FSize)
title('average epidemic')
set(gca,'Fontsize',FSize,'Xlim',[0 61]); 
box on;grid on

%INSERT YOUR CODE HERE to add plotting of other elements and optional
%things like movies