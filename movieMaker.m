%clear all data and windows
clear
clc
close all
clf

%load the data we need for plotting
load EnvironmentalForcing.mat
load simulationData.mat

for i = 1:length(tspan)
    %skip every second set of data points so we can make this process go
    %quicker
    if(mod(i,1) == 0)
        %clear the old figure
        clf
        FSize = 14; %fontsize for plots    
        %plot the lines
        plot(tspan(1:i),S_ave(1:i),'-k','LineWidth',2);
        hold on
        plot(tspan(1:i),L_ave(1:i),'-r','LineWidth',2);
        plot(tspan(1:i),I_ave(1:i),'-b','LineWidth',2);
        plot(tspan(1:i),R_ave(1:i),'-g','LineWidth',2);
        plot(tspan(1:i),P_ave(1:i),'-c','LineWidth',2);
        plot(tspan(1:i),E_ave(1:i),'-y','LineWidth',2);
        plot(tspan(1:i),B_ave(1:i),'-m','LineWidth',2);
        plot(tspan(1:i),F_ave(1:i),'-k','LineWidth',2);        
        %plot the current point
        plot(tspan(i),S_ave(i),'ko','LineWidth',3)
        plot(tspan(i),L_ave(i),'ro','LineWidth',3)
        plot(tspan(i),I_ave(i),'bo','LineWidth',3)
        plot(tspan(i),R_ave(i),'go','LineWidth',3)
        plot(tspan(i),P_ave(i),'co','LineWidth',3)
        plot(tspan(i),E_ave(i),'yo','LineWidth',3)
        plot(tspan(i),B_ave(i),'mo','LineWidth',3)
        plot(tspan(i),F_ave(i),'ko','LineWidth',3)
        %label the plot
        legend({'Susceptible','Latent','Infected','Recovered','Population','External','Berries','Spreading'},'Location','NorthWest');
        xlabel('time (days)','Fontsize',FSize);
        xlim([0 60]);
        ylabel('Population (fraction of initial)','Fontsize',FSize);
        ylim([0 1]);
        title('average epidemic')
        set(gca,'Fontsize',FSize,'Xlim',[0 61]); 
        box on;
        grid on;    
        %save the current figure to a vector
        MV(i) = getframe;
    else
        %skip every second data point
    end
end
%create a writer object
writer = VideoWriter('Population.mp4','MPEG-4');
%set frame rate
writer.FrameRate=30;
%open the writer object and write all the frames to an mp4 file
open(writer);
writeVideo(writer,MV);
close(writer);



