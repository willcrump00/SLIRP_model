function [detect_flag,diseased_plant,cost] = scouting_routine(t_idx)
global cost;
global detect_flag
global searches
%if already found we dont need to run this script and update cost
if detect_flag ==1
    return;
end

%no scouting before day 10 because very low chance of infected plants due
%to latent period - so dont run
if t_idx < 239
    return;
end
%for days 10 through 15 we run a few drones in random location
%scouting algy for days 10-15
%we search plants 240 through 360 on these days
if t_idx == 239 || t_idx == (239+24) || t_idx == (239+48) || t_idx == (239+72) || t_idx == (239+96) || t_idx == (239 + 120)
movement_speed = .00056; %m/s
%very small detect size appropriate for early infections
%max of 2 plants searched per drone
detect_size = 20 * movement_speed; %mm
%we are going to psuedo random search
for i = 1:24
    %24 is the number of plants searched by 12 drones in 1 hr at mv speed
    searched_plants(i) = t_idx + i+searches;
    infection_size = vine(t_idx+i+searches).I;
    %check if infection exists and our detect size is smaller than
    %infection size
    if infection_size > detect_size && vine(t_idx+i+searches).IsInfect == 1
        detect_flag = 1;
        diseased_plant = (t_idx+i+searches);
        break;
    end

end
%searches global adds 1 to index so that we dont search the same plant as
%we ended our search with on previous day
searches = searches +1;

%updating cost for 12x drones searching for approx 1 hour
cost = cost +  (12*100);
end

%after day 15 more t_idx  ==383 --> 503


% confirming detection success
for i = 1:length(scanned_vines)
if dectect_size > vine(scanned_vines(i)).I
    detect_flag = 1;
    diseased_plant = scanned_vines(i);
end
end


%adding cost if infected vines are not found
infected_sum = sum(vine.IsInfect);
if detect_flag == 0 && infected_sum>0
    cost = cost + 1000*infected_sum;
end

end
