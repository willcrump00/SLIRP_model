%%% Function to compute the deposition to a point (or set of points)
%%% using a version of the Gaussian plume model described by Prussian et 
%%% al. (2015 Ag&For).  As INPUTS the function REQUIRES the X and Y 
%%% positions to calculate the concentration at (as X,Y point coords in 
%%% meters from the source location), the height of the canopy (H), the 
%%% wind speed (m/s), and the wind direction (degrees).
%%%
%%% Example function call with minimum parameters:
%%%
%%% [C] = GaussianPlumeDep(X,Y,WindSpeed,WindDir,dep_area,Q)
%%%

function [C] = GaussianPlumeDep(X,Y,WindSpeed,WindDir,dep_area,Q)

%parameters for simple spread rate model (following Miller et al., 2018
%perpendicular spread fits)
sigy0 = 1.0;
my = 1.0;
set_vel=0.1; %accounts for gravity and deposition to leaves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculate the 'outgoing' wind direction in radians
WindDir=WindDir*pi/180;
dep_length=sqrt(dep_area/10000); %convert to m^2

%initialize deposition array
C = nan(1,length(X));
%calculate distance in rotated coords
x = X*cos(WindDir)+Y*sin(WindDir);
PosX = (x>0) ; %skip if its invalid (Gauss Plumes are invalid in the -dir)
y = -X*sin(WindDir)+Y*cos(WindDir);
sigy = sqrt(sigy0^2+my^2*(x(PosX)).^2);

C(PosX) = (1./(2*3.14*sigy*WindSpeed)).*exp(-((y(PosX)).^2)./(2*sigy.^2));
C(PosX) = exp(-0.05*abs(x(PosX))).*C(PosX); %Miller et al (2018) removal rate
C(PosX) = real(set_vel.*C(PosX)*dep_length.*Q);
%note that this C is the deposition in amount/second (a rate)
end