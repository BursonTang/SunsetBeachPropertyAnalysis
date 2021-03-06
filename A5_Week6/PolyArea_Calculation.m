%% Description

%{
A numerical method is used to calculate the Beach volume, which is
compared with the Surface Volume Tool in ArcGIS
%}

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
packx = NaN(24,25,2500);   % Pack for distance from -8isobath
packz = NaN(24,25,2500);   % Pack for corresponding elevation
packh = zeros(24,25,2500);   % Pack for corresponding height from MSL
packt = zeros(24,25);      % Pack for counting the number of data pairs

%% Convert data to data pack

% Initialize counting variables for dem, line and data points
nd = 1; nl = 1; np = 1; t = length(rawdata(:,3));

for a = 1: t
    
    % Convert data from raw data to data pack
    packx(nd,nl,np) = rawdata(a,1);
    packh(nd,nl,np) = rawdata(a,2)-0.798;
    
    % For z value
    if rawdata(a,2) == -99 && rawdata(a,1)<=800
        packz(nd,nl,1:np) = NaN;
    end
    
    if rawdata(a,2) < 0 %|| rawdata(a,1) < 700
        packz(nd,nl,np) = NaN; % Set NaN
    else
        packz(nd,nl,np) = rawdata(a,2);
    end
    
    % For h value
    if packh(nd,nl,np)<0
        packh(nd,nl,np) = 0;
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

%% Eliminate Wave Data
% for c = 1:20   % 24 dem datasets
%     
%     for n = 1 : 25   % 25 crosssections
%         for i = packt(c,n):-1:2
%             if ( packz(c,n,i-1) > packz(c,n,i) && packz(c,n,i)<0.6 )% || ( isnan(packz(c,n,i)) && packz(c,n,i) < 3 )
%                 packz(c,n,1:i)= NaN;
%                 break
%             end
%         end
%     end
%     
% end

%% Visualization Part 1

% Plot graphs on all valid dem dataset (20 in total) at crosssection#n
% for n = 1:25;
n = 25;

% Variable recording data from each dem file for debugging use
x = NaN(2500,24);
z = NaN(2500,24);

figure(1);grid on; hold on;

for b = 1: 15
    x(1:packt(b,n),b) = squeeze(packx(b,n,1:packt(b,n)));
    z(1:packt(b,n),b) = squeeze(packz(b,n,1:packt(b,n)));
    if b<= 5
        plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
            ,'-','LineWidth',2)
    else if b<=10
            plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                ,':','LineWidth',2)
        else if b<=15
                plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                    ,'-.','LineWidth',2)
            else
                plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                    ,'--','LineWidth',2)
            end
        end
    end
end

% Extract the maximum distance from the farthest point of the longest
% profile to the -8 isobath

leftb = [720,750,800,850,900,900,850,850,820,850,850,850,750,920,920,920,930,980,980,950,920,900,900,880,850];
xlim([leftb(n),nanmax(nanmax(squeeze(packx(:,n,:))))+10])

title(['Profiles at Cross Section ',int2str(n)])
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'...
    ,'Location','SouthEast')
xlabel('Distance from -8 Isobath (m)')
ylabel('Elevation (m)')

% % % Debug part
% pause (10)
% close all
% end

%% Visualization Part 2
% Extract single profile graph on dem#m at crosssection#n
figure(2);grid on; hold on;
m = 14;   % DEM#m, the x and z data for the single line is x or z(m,1:2500)

plot(squeeze(packx(m,n,1:packt(m,n))),squeeze(packz(m,n,1:packt(m,n))))
xlim([leftb(n),nanmax(x(:,b)+10)])

title(['Profiles at Cross Section ',int2str(n),' on the DEM dataset ',int2str(m)])
xlabel('Distance from -8 Isobath (m)')
ylabel('Elevation (m)')

%% Calculate the Area Under Curve (Above MSL)
Area = zeros(15,25);

for nd = 1:15
    for nl = 1:25
        for np = 2: packt(nd,nl)
            height_m = (packh(nd,nl,np-1) + packh(nd,nl,np))/2;
            dis_x = packx(nd,nl,np) - packx(nd,nl,np-1);
            Area_sub = dis_x*height_m;
            Area(nd,nl) = Area(nd,nl) + Area_sub;
        end
    end
end

dis_alongb = [95,96,106,102,91,110,102,84,103,108,101,98,92,106,92,106,97,...
    104,111,101,85,114,103,90];

Volume = zeros(1,15);

for nd = 1:15
    for nl = 2:25
        Area_m = (Area(nd,nl-1) + Area(nd,nl))/2;
        Volume(nd) = Volume(nd) + Area_m*dis_alongb(nl-1);
    end
end