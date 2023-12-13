function [detect_flag,cost,infected_plant,latent_plant] = scouting_routine(tspan,t_idx,vine)
global detect_flag cost drone_last infected_plant p q latent_plant %my added globals for scouting function
A = 5000;
%movement vars each is 1 step in respective direction
L = -1;
R = 1;
U = 50;
D = -50;

%if already found we dont need to run this script and update cost
%if detect_flag ==1
    %fprintf('\n diseased plant found number: %.2f \n',infected_plant);
   % return;
%end

%no scouting before day 13 because very low chance of infected plants due
%to latent period - its actually cheaper to let any early infections fester
%so we can more accurately detect
%starting on day 13 search all plants daily
%2 hours before the end of the day to maximize in fection detection rates
if t_idx < 241
    cost = 0;
    detect_flag = 0;
    drone_last = [1225,1226,1276,1275,1274,1224,1174,1176,1,2451,2500,50];
    p = 0;
    q=0;
    latent_flag=0;
    %giving default values from within the function because it didnt like
    %assigning them fr script
end
%days 12-16 searching approx 84 plants x day
if t_idx >=311  && t_idx <= 383 && mod((round(tspan(t_idx),4) - 0.9167),1) == 0
    if detect_flag ==1
        return
    end
    if t_idx == 383 %reverse directions & add offset
        D=-D;
        R=-R;
        U=-U;
        L=-L;
    end
    %min infection size is about 1.34*10^-5 cm^2
    cost = cost +1200; %cost of renting 12 drones for 1 hr
    movement_speed = .00195; %m/s = 7 plants per drone
    detect_size = (20 * movement_speed)/(2*10); %mm (convert to cm & radius)
    detect_area = (detect_size^2)*pi; %cm ^2   
%spiral pattern reaching out from center with 4 drones in each corner as
%well
for j = 1:12
    searched_plants(j,1) = drone_last(j);%technically we dont have to move at all to search the first plants
    for i = 2:8
        if j == 1
            searched_plants(j,i) = searched_plants(j,i-1) + D;
        end
        if j == 2
            searched_plants(j,i) = searched_plants(j,i-1) + R;
        end
        if j ==3 
             searched_plants(j,i) = searched_plants(j,i-1) + R+U;
        end
        if j ==4
             searched_plants(j,i) = searched_plants(j,i-1) + U;
        end
        if j ==5
             searched_plants(j,i) = searched_plants(j,i-1) + L+U;
        end
        if j ==6
             searched_plants(j,i) = searched_plants(j,i-1) + L;
        end
        if j == 7
             searched_plants(j,i) = searched_plants(j,i-1) + L + D;
        end
        if j == 8
             searched_plants(j,i) = searched_plants(j,i-1) + R+D;
        end
        if j ==9
             searched_plants(j,i) = searched_plants(j,i-1) + U;
        end
        if j ==10
             searched_plants(j,i) = searched_plants(j,i-1) + R;
        end
        if j ==11
             searched_plants(j,i) = searched_plants(j,i-1) + D;
        end
        if j ==12
             searched_plants(j,i) = searched_plants(j,i-1) + L;
        end

        %testing plants as we go
        if(searched_plants(j,i)>2500)
            fprintf('\n out of bounds at %.1f',searched_plants(j,i))
            error('\n we exceeded bounds at %.1f j and %.1f i',j,i)
        end
        infection_size = vine(searched_plants(j,i)).I(t_idx) * A; %convert infection size of each plant to area of infection on plant
        %check if infection exists and our detect size is smaller than
        %infection size
        if infection_size ~= 0
            fprintf('\n infection area: %.6f',infection_size)
            fprintf('\n detect area: %.6f',detect_area)
            if infection_size>detect_area
                latent_flag =1;
                q=q+1;
                latent_plant(q) =searched_plants(j,i);
            end

            if infection_size > detect_area && vine(searched_plants(j,i)).IsInfect == 1
                detect_flag = 1;
                p=p+1;
                infected_plant(p) = searched_plants(j,i);
                fprintf('\n diseased plant found, locating and exterminating plant \n number: %.2f',infected_plant(p));
                
            break;
            end

        end

    end
    drone_last(j) = searched_plants(j,i); % ensure this is implemented correctly (good)
    %sets pos of each drone for last day
end
end


%day 17 search fast and gridded
if t_idx == 407 || t_idx == 455
    if detect_flag ==1
        return
    end
    %min if theres one thats infected for 3 days by now  - size is about
    %,055cm^2
    cost = cost +1200; %cost of renting 12 drones for 1 hr
    movement_speed = .0423; %m/s = 152 plants per drone
    detect_size = (20 * movement_speed)/(2*10); %mm (convert to cm & radius)
    detect_area = (detect_size^2)*pi; %cm ^2
    %new drone last because new pattern
    drone_last = [1,151,301,451,601,751,901,1051,1201,1351,1451,1601];
    %this pattern should be easy , go across then up 1 then back then up 1
    %then across = 152 plants
   for j = 1:12
    searched_plants(j,1) = drone_last(j); %technically we dont have to move at all to search the first plants
    i=2; % reset indexing marker
        while (i<51)
        searched_plants(j,i) = searched_plants(j,i-1)+R;
        i=i+1;
        end
        searched_plants(j,51) = searched_plants(j,i-1)+U;
        i=i+1;
        while (i<101)
        searched_plants(j,i) = searched_plants(j,i-1)+L;
        i=i+1;
        end
        searched_plants(j,101) = searched_plants(j,i-1)+U;
        i = i+1;
        while (i<152)
        searched_plants(j,i) = searched_plants(j,i-1)+R;
        i=i+1;
        end
        %testing plants as we go
        for i = 1:length(searched_plants)
        infection_size = vine(searched_plants(j,i)).I(t_idx) * A; %convert infection size of each plant to area of infection on plant
        %check if infection exists and our detect size is smaller than
        %infection size
        if infection_size ~= 0
            fprintf('\n infection area: %.6f',infection_size)
            fprintf('\n detect area: %.6f',detect_area)
             if infection_size>detect_area
                latent_flag =1;
                q=q+1;
                latent_plant(q) =searched_plants(j,i);
            end
            if infection_size > detect_area && vine(searched_plants(j,i)).IsInfect == 1
                detect_flag = 1;
                p = p+1;
                infected_plant(p) = searched_plants(j,i);
                fprintf('\n diseased plant found, locating and exterminating plant \n number: %.2f',infected_plant(p));
            break;
            end

        end
        end
   end

%day 18 same speed but starting from top or day 20 if none found
if t_idx == 431 || t_idx == 479
    if detect_flag ==1
        return
    end
    %min if theres one thats infected for 3 days by now  - size is about
    %,055cm^2
    cost = cost +1200; %cost of renting 12 drones for 1 hr
    movement_speed = .0423; %m/s = 152 plants per drone
    detect_size = (20 * movement_speed)/(2*10); %mm (convert to cm & radius)
    detect_area = (detect_size^2)*pi; %cm ^2
    %new drone last because new pattern
    drone_last = [2451,2301,2151,2001,1851,1701,1501,1351,1201,1051,901,751];
    %this pattern should be easy , go across then down then across
    %then across = 152 plants

   for j = 1:12
    searched_plants(j,1) = drone_last(j); %technically we dont have to move at all to search the first plants
    i=2; % reset indexing marker
        while (i<51)
        searched_plants(j,i) = searched_plants(j,i-1)+R;
        i=i+1;
        end
        searched_plants(j,51) = searched_plants(j,50)+D;
        i=i+1;
        while (i<102)
        searched_plants(j,i) = searched_plants(j,i-1)+L;
        i=i+1;
        end
        searched_plants(j,102) = searched_plants(j,101)+D;
        i = i+1;
        while (i<153)
        searched_plants(j,i) = searched_plants(j,i-1)+R;
        i=i+1;
        end
        %testing plants as we go
        for i = 1:length(searched_plants)
        infection_size = vine(searched_plants(j,i)).I(t_idx) * A; %convert infection size of each plant to area of infection on plant
        %check if infection exists and our detect size is smaller than
        %infection size
        if infection_size ~= 0
            fprintf('\n infection area: %.6f',infection_size)
            fprintf('\n detect area: %.6f',detect_area)
             if infection_size>detect_area
                latent_flag =1;
                q=q+1;
                latent_plant(q) =searched_plants(j,i);
            end
            if infection_size > detect_area && vine(searched_plants(j,i)).IsInfect == 1
                detect_flag = 1;
                 p = p+1;
                infected_plant(p) = searched_plants(j,i);
                fprintf('\n diseased plant found, locating and exterminating plant \n number: %.2f',infected_plant(p));
            break;
            end

        end
        end
   end
end

%after day 18 if nothing we wait because its probably a late starting
%infection - then repeat the day 17 and 18 procedures for day 19-20

%day 20 and no infection!!! oh shit



%adding cost if infected vines are not found

infected_sum = sum([vine.IsInfect]);
if detect_flag == 0 && infected_sum>0 && t_idx >241 && mod(tspan(t_idx),1)==0
    cost = cost + 1000*infected_sum;
end
infected_plant;
latent_plant;
end
