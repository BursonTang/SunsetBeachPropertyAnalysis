%% Description
% The code is for a CEE 199 assignment, this code is pretty similart to the
% script "Beach_summit_slope_v1" at A4 file. Only wave eliminated part
% is modified

% Also, this Script add a visulization part for slope and crest elevation.
% The mean value in along longshore direction is plotted

% 24 DEM are used in arcgis to extrac cross shore profile of the Sunset
% Beach. The last 4 DEM data are considered as invalid and thus only 20 DEM
% are used to plot

% The code eliminate the data affected by wave height, and the maximum
% foreshore elevation and foreshore beach slope are caculated in this code

clear all; close all; clc

%% Initialization
% Read CSV files
rawdata = csvread('LinePack_modified_2.csv',1,1);

% Mean sea level and mean higher high water (Los Angeles Tide Gauge DY1083)
navd = .063;
msl = .861 - navd;
mhhw= 1.674 - navd;

% Preallocate for data, 24 DEMs all have 25 lines, which contain less than
% 1000 data points
dem = 24; line = 25;
packx = NaN(24,25,2500);   % Pack for distance from -8isobath
packz = NaN(24,25,2500);   % Pack for corresponding elevation
packt = zeros(24,25);      % Pack for counting the number of data pairs

%% Convert data to data pack

% Initialize counting variables for dem, line and data points
nd = 1; nl = 1; np = 1; t = length(rawdata(:,3));

for a = 1: t
    
    % Convert data from raw data to data pack
    packx(nd,nl,np) = rawdata(a,1);
    
    if rawdata(a,2) == -99 && rawdata(a,1)<=800
        packz(nd,nl,1:np) = NaN;
    end
    
    if rawdata(a,2) < 0 || rawdata(a,1) < 700
        packz(nd,nl,np) = NaN; % Set NaN
    else
        packz(nd,nl,np) = rawdata(a,2);
    end
    np = np + 1;
    
    if a < t
        if rawdata(a+1,3) - rawdata(a,3) == 1
            packt(nd,nl) = np-1;
            nl = nl+1;
            np = 1;
            
        else if rawdata(a+1,3) - rawdata(a,3) ~= 0
                packt(nd,nl) = np-1;
                nd = nd+1;
                nl = 1;
                np = 1;
                
            end
        end
    end
end

% Eliminate Wave Data
for c = 1:20   % 24 dem datasets
    
    for n = 1 : 25   % 25 crosssections
        for i = packt(c,n):-1:2
            if ( packz(c,n,i-1) > packz(c,n,i) && packz(c,n,i)<0.798 )% || ( isnan(packz(c,n,i)) && packz(c,n,i) < 3 )
                packz(c,n,1:i-1)= NaN;
                break
            end
        end
    end
    
end


%% Find the maximum beach elevation and calculate the beach slope
max_BE = -99*ones(24,25);

for n = 1:25   % 25 crosssections
    
    for c = 1 : 20   % 24 dem datasets
        % Start from points with z that is close to 2.5
        [~,i]= nanmin(abs(squeeze(packz(c,n,:))-2.5));
        
        while max_BE(c,n) < packz(c,n,i)
            max_BE(c,n) = packz(c,n,i);
            i = i + 1;
        end
    end
end

%% Calculate the Beach Slope
% Preallocate for start point and end point for the slope calculation
s_b = ones(24,25);
e_b = ones(24,25);
count = 0;
% find the start point
for n = 1 : 25   % 25 crosssections
    for c = 1 : 20   % 24 dem datasets

        % Use the highest elevation points as the end points
        midn = find ( packz(c,n,:) == max_BE(c,n) );
        e_b(c,n) = midn(1);
        
        
        for i = 2 : e_b(c,n)-1
            if packz(c,n,i+1) <= packz(c,n,i) || isnan(packz(c,n,i))
            s_b(c,n) = i+1;
            if isnan(packz(c,n,i+1))
                s_b(c,n) = i+2;
            end
            end
        end

    end
end

% Calculation of beach slope
slope_b = NaN(24,25);
for n = 1:25   % 25 crosssections
    for c = 1 : 15   % 24 dem datasets
        % Get the elevation and distance for both points
        z1 = packz(c,n,s_b(c,n));
        x1 = packx(c,n,s_b(c,n));
        
        z2 = packz(c,n,e_b(c,n));
        x2 = packx(c,n,e_b(c,n));
        
        slope_b(c,n) = (z2-z1)/(x2-x1);

    end
end

%% Visualization Part 

dis_x = [95,96,106,102,91,110,102,84,103,108,101,98,92,106,92,106,97,...
    104,111,101,85,114,103,90];
x_ab = zeros(1,25);

% Summit Elevation
figure(1);grid on; hold on;

for nl = 2:25 
    x_ab(nl) = x_ab(nl-1) + dis_x(nl-1);
end

% for nd = 1:15
%     if nd<=5
%     plot(x_ab,max_BE(nd,:),'-','LineWidth',1.5);
%     else if nd<=10
%             plot(x_ab,max_BE(nd,:),'--','LineWidth',1.5);
%         else plot(x_ab,max_BE(nd,:),':','LineWidth',1.5);
%         end
%     end
% end

plot(x_ab,nanmean(max_BE(1:15,:)))

xlabel('Along Beach Distance(m)')
ylabel('Foreshore Maximum Elevation(m)')
title('Foreshore Maximum Elevation')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'...
    ,'Location','NorthEast')
xlim([-10,2410])
%set (gca,'Ydir','reverse')

% MHHW shoreline
figure(2);grid on; hold on;

% for nd = 1:15
%     if nd<=5
%     plot(x_ab,slope_b(nd,:),'-','LineWidth',1.5);
%     else if nd<=10
%             plot(x_ab,slope_b(nd,:),'--','LineWidth',1.5);
%         else plot(x_ab,slope_b(nd,:),':','LineWidth',1.5);
%         end
%     end
% end

plot(x_ab,nanmean(slope_b(1:15,:)))

xlabel('Along Beach Distance(m)')
ylabel('Slope')
title('Slope on Different Datasets at Different Transects')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'...
    ,'Location','NorthEast')
xlim([-10,2410])
%set (gca,'Ydir','reverse')