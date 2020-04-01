%% Description
% The script is one of the final scripts CEE 199. 
% Error bar is added for the RA assignment
clear all; close all; clc

%% Initialization
% Read CSV files
rawdata = csvread('TransectsData_afex_modified.csv',1,1);

% Mean sea level and mean higher high water (Los Angeles Tide Gauge DY1083)
navd = .063;
msl = .861 - navd;
mhhw= 1.674 - navd;

% Preallocate for data, 17 DEMs all have 25 lines, which contain less than
% 1000 data points
packx = NaN(17,25,2500);   % Pack for distance from -8isobath
packz = NaN(17,25,2500);   % Pack for corresponding elevation
packt = zeros(17,25);      % Pack for counting the number of data pairs

%% Convert data to data pack

% Initialize counting variables for dem, line and data points
nd = 1; nl = 1; np = 1; t = length(rawdata(:,3));

for a = 1: t
    
    % Convert data from raw data to data pack
    packx(nd,nl,np) = rawdata(a,1);
    
    if rawdata(a,2) == -99 && rawdata(a,1)<=800
        packz(nd,nl,1:np) = NaN;
    end
    
    if rawdata(a,2) < 0 || rawdata(a,1) < 700 || rawdata(a,2) > 10
        packz(nd,nl,np) = NaN; % Set NaN
    else
        packz(nd,nl,np) = rawdata(a,2);
    end
    
    % Fix unvalid high elevation point
    if rawdata(a,2) > 10
        packz(nd,nl,np) = (rawdata(a-1, 2) + rawdata(a+1, 2))/2;
    end
    
    np = np + 1;
    
    if a < t
        if rawdata(a+1,3) - rawdata(a,3) == 1
            packt(nd,nl) = np-1;
            nl = nl+1;
            np = 1;
            
        else if rawdata(a+1,3) - rawdata(a,3) == -24
                packt(nd,nl) = np-1;
                nd = nd+1;
                nl = 1;
                np = 1;
                
            end
        end
    else
        packt(nd,nl) = np-1;
    end
end

% Eliminate Wave Data
for c = 1:17   % 17 dem datasets
    
    for n = 1 : 25   % 25 crosssections
        for i = packt(c,n):-1:2
            if ( packz(c,n,i-1) >= packz(c,n,i) && packz(c,n,i)<1.40 )% || ( isnan(packz(c,n,i)) && packz(c,n,i) < 3 )
                %1.6110 or 0.798?  Set as 1.40 to include MHHW position
                packz(c,n,1:i-1)= NaN;
                break
            end
        end
    end
    
end

%% Visualize the profile of transects

% Plot graphs on all valid dem dataset (20 in total) at crosssection#n
% for n = 1:25;
n = 23;

% Variable recording data from each dem file for debugging use
x = NaN(2500,17);
z = NaN(2500,17);

figure(1);grid on; hold on;

for b = 1: 17
    x(1:packt(b,n),b) = squeeze(packx(b,n,1:packt(b,n)));
    z(1:packt(b,n),b) = squeeze(packz(b,n,1:packt(b,n)));
    if b<= 4
        plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
            ,'-','LineWidth',2)
    else if b<=8
            plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                ,':','LineWidth',2)
        else if b<=12
                plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                    ,'--','LineWidth',2)
            else
                plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
                    ,'-.','LineWidth',2)
            end
        end
    end
end

% Extract the maximum distance from the farthest point of the longest
% profile to the -8 isobath

leftb = fliplr([720,750,800,850,900,900,850,850,820,850,850,850,750,...
    920,920,920,930,980,980,950,920,900,900,880,850]);
xlim([leftb(n),nanmax(nanmax(squeeze(packx(:,n,:))))+10])

title(['Beach Profiles at Transect ',int2str(n)])

set(gca,'FontSize',20)
legend({'Sep 1997','Apr 1998','Sep 2004','Apr 2005','Oct 2005','Mar 2006',...
    'Oct 2006','Apr 2007','Dec 2007','Apr 2008','Oct 2008','Mar 2009',...
    '2009~2011','Sep 2014','Mar 2016','May 2016','Oct 2016'},...
    'Location','SouthEast','Fontsize',10)

xlabel('Cross Shore Distance from 8-meter Isobath (m)')
ylabel('Elevation NAVD88 (m)')
% f1 = gca;

% % Debug part
% % pause (5)
% % close all
% % end

%% Find the maximum foreshore elevation and calculate the beach slope

% Maximum foreshore elevation
max_BE = -99*ones(17,25);

% Maximum beach elevation
max_BE_2 = nan(17,25);

for n = 1:25   % 25 crosssections
    
    for c = 1 : 17   % 17 dem datasets
        
        % Start from points with z that is close to 2.5
        [~,i]= nanmin(abs(squeeze(packz(c,n,:))-2.5));
        
        while max_BE(c,n) < packz(c,n,i)
            max_BE(c,n) = packz(c,n,i);
            i = i + 1;
        end
        
        % Find the maximum beach elevation (Odelette's method)
        max_BE_2(c,n) = nanmax ( squeeze( packz(c,n,:) ) );
    end
end

%% Calculate the Beach Slope

% Preallocate for start point and end point for the slope calculation
s_b = ones(17,25);
e_b = ones(17,25);


% find the starting and ending point
for n = 1 : 25   % 25 crosssections
    for c = 1 : 17   % 17 dem datasets

        % Use the highest elevation points as the end points
        midn = find ( packz(c,n,:) == max_BE(c,n) );
        e_b(c,n) = midn(1);   % Use the first crest
        
        % Find the start point of beach
        for i = 2 : e_b(c,n)-1
            if  packz(c,n,i+1) <= packz(c,n,i) || isnan(packz(c,n,i)) ||...
                    isnan(packz(c,n,i+1)) 
            s_b(c,n) = i+1;
%             if isnan(packz(c,n,i+1))
%                 s_b(c,n) = i+2;
%             end
            end
        end

    end
end

% Calculation of beach slope
slope_b = NaN(17,25);
for n = 1:25   % 25 crosssections
    for c = 1 : 17   % 17 dem datasets
        % Get the elevation and distance for both points
        z1 = packz(c,n,s_b(c,n));
        x1 = packx(c,n,s_b(c,n));
        
        z2 = packz(c,n,e_b(c,n));
        x2 = packx(c,n,e_b(c,n));
        
        slope_b(c,n) = (z2-z1)/(x2-x1);

    end
end

%% Calculate mean foreshore slope and mean maximum crest elevation

% Temporally averaged slope and maximum crest elevation
mbc = nanmean(max_BE);
mbs = nanmean(slope_b); % slope~ 1:mbs ~ x:z

% Mean slope and maximum crest among 17*25 value
t_avg_bs = nanmean(mbs);
t_avg_bc = nanmean(nanmean(max_BE_2));

%% Error Bar

erb_maxBE = max(max_BE) - mbc;
erb_minBE = min(max_BE) - mbc;

erb_maxBS = max(slope_b) - mbs;
erb_minBS = min(slope_b) - mbs;
% erb_maxBE = max(max_BE);
% erb_minBE = min(max_BE);

%% Plot Maximum Crest Elevation
figure(2);grid on; hold on;

x_ab = 1:25;

yyaxis left
plot(x_ab,mbc,'-r','LineWidth',1.5)
errorbar(x_ab,mbc,erb_minBE,erb_maxBE,'r')

% Plot the maximum beach elevation
% plot(x_ab,nanmean(max_BE_2),'--','LineWidth',1.5)

ylabel('Foreshore Maximum Elevation(ft NAVD88)')
% xlabel('Along Beach Distance(m)')
% title('Mean Foreshore Maximum Elevation')
% legend('Maximum Elevation','Location','NorthEast')
% xlim([-10,2410])
% set (gca,'Ydir','reverse')

%% Plot Foreshore Slope

yyaxis right
plot(x_ab,mbs,'-b','LineWidth',1.5)
errorbar(x_ab,mbs,erb_minBS,erb_maxBS,'b')

xlabel('Transect Number')
ylabel('Slope')
title('Mean Slope and Foreshore Maximum Elevation at Different Transects')
legend({'Elevation','The Error Bar of Elevation','Slope'},'Location','NorthEast')
set(gca,'FontSize',20)
% legend('Elevation','e','Slope','e','Location','NorthEast')
% xlim([-50,2450])
xlim([0,26])
set(gca,'XTick',1:25)

%% Visualization Part 2, time

% % The collceted date of new 17 datasets(week7),
% Datestring = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
%     '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';...
%     '09/2014';'03/2016';'05/2016';'10/2016'];
% 
% formatin = 'mm/yyyy';
% ta = datenum(Datestring,formatin);
% 
% mtc = nanmean(max_BE,2);
% mts = nanmean(slope_b,2);
% 
% % % Error Bar
% % maxc = nanmax(max_BE);
% % minc = nanmin(max_BE);
% % 
% % maxs = nanmax(slope_b);
% % mins = nanmin(slope_b);
% %% Summit Elevation
% figure(2);grid on; hold on;
% 
% yyaxis left
% plot(ta,[mtc(1:12);mtc(14:17)])
% % errorbar(x_ab,nanmean(max_BE),maxc,minc)
% 
% ylabel('Foreshore Maximum Elevation(m)')
% % xlabel('Along Beach Distance(m)')
% % title('Mean Foreshore Maximum Elevation')
% % legend('Maximum Elevation','Location','NorthEast')
% % xlim([-10,2410])
% % set (gca,'Ydir','reverse')
% 
% %% Foreshore Elevation
% 
% yyaxis right
% plot(ta,[mts(1:12);mts(14:17)])
% % errorbar(x_ab,nanmean(slope_b),maxs,mins)
% 
% dateFormat = 'mm/yyyy';
% datetick('x',dateFormat)
% 
% xlabel('Time')
% ylabel('Slope')
% title('Mean Slope and Foreshore Maximum Elevation at Different Transects')
% legend('Elevation','Slope','Location','NorthEast')
% % legend('Elevation','e','Slope','e','Location','NorthEast')
% %xlim([-50,2450])
