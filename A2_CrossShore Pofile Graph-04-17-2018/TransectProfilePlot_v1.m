%% Description

%{ 
28 cross-shore transects profiles are firsted extracted from three Digital
Elevation Model (DEM) data with ArcGIS. The profile data is stored in the
csv files, in which first column represents the cross-shore distance
starting from the onshore edges, and second column indicating the
elevation.

Question: Is it possible to make the distance accumulation starts from a
line, which can be assumed to be stationary over time?
=> find the -8m isobath line, and make the distance of profile lines
starting from the -8m isobath line.
%}
clear all;close all;clc;

%%  Read Files and Preallocate for
% Read CSV files
data2009 = csvread('2009.csv',1,1);
data2014 = csvread('2014.csv',1,1);
data2016 = csvread('2016.csv',1,1);

% Count profiles number
pronum = length(find(data2009(:,1)==0));

% Preallocate matrix for

% Elevation data
y09 = NaN(pronum,1000);
y14 = NaN(pronum,1000);
y16 = NaN(pronum,1000);

% Distance data
x09 = NaN(pronum,1000);
x14 = NaN(pronum,1000);
x16 = NaN(pronum,1000);

% The length of cross-shore line
inter09 = ones(28,1);
inter14 = ones(28,1);
inter16 = ones(28,1);

%% The loop for data of 2009
x = 1; y = 2; y09(1,1) = data2009(1,2); x09(1,1) = 0;
for i09 = 2: length(data2009(:,1))
    
    if data2009(i09,1) ~= 0
        y09(x,y) = data2009(i09,2);
        x09(x,y) = data2009(i09,1);
        inter09(x,1) = inter09(x,1)+1;
        y = y + 1;
    else % For a new profile, reset the default values
        x = x +1;
        y09(x,1) = data2009(i09,2);
        x09(x,1) = 0;
        y = 2;
    end
    
end

%% The loop for data of 2014
x = 1; y = 2; y14(1,1) = data2014(1,2); x14(1,1) = 0;
for i14 = 2: length(data2014(:,1))
    
    if data2014(i14,1) ~= 0
        y14(x,y) = data2014(i14,2);
        x14(x,y) = data2014(i14,1);
        inter14(x,1) = inter14(x,1)+1;
        y = y + 1;
    else % For a new profile, reset the default values
        x = x +1;
        y14(x,1) = data2014(i14,2);
        x14(x,1) = 0;
        y = 2;
    end
    
end

%% The loop for data of 2016
x = 1; y = 2; y16(1,1) = data2016(1,2); x16(1,1) = 0;
for i16 = 2: length(data2016(:,1))
    
    if data2016(i16,1) ~= 0
        y16(x,y) = data2016(i16,2);
        x16(x,y) = data2016(i16,1);
        inter16(x,1) = inter16(x,1)+1;
        y = y + 1;
    else % For a new profile, reset the default values
        x = x +1;
        y16(x,1) = data2016(i16,2);
        x16(x,1) = 0;
        y = 2;
    end
    
end

% plot the reference map
figure(1)
Sunset = imread('Sunset Beach.png');
imshow(Sunset);

% Prompt to get the desired Line Order
LO = input('While profile graph(0~27) do you want to plot?');
LO = LO+1;


% Plot
figure(2); hold on; grid on;
plot(x09(LO,1:inter09(LO)),y09(LO,1:inter09(LO)),'r','LineWidth',2)
plot(x14(LO,1:inter14(LO)),y14(LO,1:inter14(LO)),'b','LineWidth',2)
plot(x16(LO,1:inter16(LO)),y16(LO,1:inter16(LO)),'g','LineWidth',2)

title('Topobathy Levels')
xlabel('Cross-shore Distance (meters)')
ylabel('Elevation (NAVD88 meters)')
legend('2009 (March) SIO Lidar: Southern CA',...
    '2014 USACE NCMP Topobathy Lidar DEM',...
    '2016 USGS West Coast El-Nino Lidar DEM','Location','SouthOutside')
%% Subplot

figure(3);hold on;
for i = 1:28
subplot(4,7,i);hold on;

plot(x09(i,1:inter09(i)),y09(i,1:inter09(i)),'r','LineWidth',2)
plot(x14(i,1:inter14(i)),y14(i,1:inter14(i)),'b','LineWidth',2)
plot(x16(i,1:inter16(i)),y16(i,1:inter16(i)),'g','LineWidth',2)
xlim([0,150])

end

% title('Topobathy Levels')
% xlabel('Cross-shore Distance (meters)')
% ylabel('Elevation (NAVD88 meters)')
% legend('2009 (March) SIO Lidar: Southern CA',...
%     '2014 USACE NCMP Topobathy Lidar DEM',...
%     '2016 USGS West Coast El-Nino Lidar DEM','Location','SouthOutside')