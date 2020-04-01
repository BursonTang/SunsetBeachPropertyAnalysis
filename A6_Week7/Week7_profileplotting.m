%% Description
% The script is one of the final scripts CEE 199. It's the debugged and 
% modified version of Week6_profileplotting.

% Invalid DEMs are excluded in this analysis. Only 17 DEMs are used.

% Cross section profiles are plotted and the graphs are modified.

% Shoreline positions are de-meaned by averaged shoreline 
% position at each dataset. Averaged demeaned shoreline position is plotted
% with error bar

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
    
    if rawdata(a,2) < 0 || rawdata(a,1) < 700 %|| rawdata(a,2) > 10 
        packz(nd,nl,np) = NaN; % Set NaN
    else
        packz(nd,nl,np) = rawdata(a,2);
    end
    
    % Fix unvalid high elevation point
%     if rawdata(a,2) > 10 
%         packz(nd,nl,np) = (rawdata(a-1, 2) + rawdata(a+1, 2))/2;
%     end
    
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

% Eliminate Wave Data and unvalid high elevation
for c = 1:17   % 17 dem datasets
    
    for n = 1 : 25   % 25 crosssections
        for i = packt(c,n):-1:2
            % Wave Data
            if ( packz(c,n,i-1) >= packz(c,n,i) && packz(c,n,i)<1.40 ) % || ( isnan(packz(c,n,i)) && packz(c,n,i) < 3 )
                %1.6110 or 0.798?  Set as 1.40 to include MHHW position
                packz(c,n,1:i-1)= NaN;
                break
            end
            
%             Invalid High Elevation
            if packz(c,n,i)>10 || (packz(c,n,i)>6 && (n ==10 || n == 17))
                packz(c,n,i) = (packz(c,n,i-1) + packz(c,n,i+1))/2;
            end
        end
    end
    
end

%% Calculate MSL and MHHW distance

% Preallocate for recording MSL and MHHW distance
x_msl = NaN(17,25);
x_mhhw = NaN(17,25);
z_d_msl = NaN(17,25); z_index_msl = NaN(17,25);
z_d_mhhw = NaN(17,25); z_index_mhhw = NaN(17,25);

for n = 1:25   % 25 crosssections
    
    % Preallocate for x and z (Set an oversized matrix, because different dem
    % have different data set size.
    x = NaN(2500,17);
    z = NaN(2500,17);
    
    for c = 1:17   % 17 dem datasets
        x(1:packt(c,n),c) = squeeze(packx(c,n,1:packt(c,n)));
        z(1:packt(c,n),c) = squeeze(packz(c,n,1:packt(c,n)));
        
        
        % find index for msl and mhhw
%         [z_d_msl(c,n), z_index_msl(c,n)] = nanmin(fliplr(abs(squeeze(packz(c,n,1:packt(c,n))) - msl)));
%         [z_d_mhhw(c,n), z_index_mhhw(c,n)] = nanmin(fliplr(abs(squeeze(packz(c,n,1:packt(c,n))) - mhhw)));
        
        [z_d_msl(c,n), z_index_msl(c,n)] = nanmin(flipud(abs(z(:,c) - msl)));
        [z_d_mhhw(c,n), z_index_mhhw(c,n)] = nanmin(flipud(abs(z(:,c) - mhhw)));
        
        % Find corresponding distance(Linear Interpolation)
        % For msl
        if z(2501-z_index_msl(c,n),c)>= msl && z(2501-z_index_msl(c,n)-1,c)< msl
            x_msl(c,n) = (msl-z(2501-z_index_msl(c,n)-1,c))*...
                (x(2501-z_index_msl(c,n),c) -  x(2501-z_index_msl(c,n)-1,c))/...
                (z(2501-z_index_msl(c,n),c) - z(2501-z_index_msl(c,n)-1,c))...
                +x(2501-z_index_msl(c,n)-1,c);
            
        else if z(2501-z_index_msl(c,n),c)<= msl && z(2501-z_index_msl(c,n)+1,c)> msl
                x_msl(c,n) = (msl-z(2501-z_index_msl(c,n),c))*...
                    (x(2501-z_index_msl(c,n)+1,c) -  x(2501-z_index_msl(c,n),c))/...
                    (z(2501-z_index_msl(c,n)+1,c) - z(2501-z_index_msl(c,n),c))...
                    +x(2501-z_index_msl(c,n),c);
                
                % The cross section profile doesn't cross the msl or mhhw
            else %if abs(z_d_msl(c,n))< 0.1  % Difference Tolerance < 0.1 meter
                    %x_msl(c,n) = x(2501-z_index_msl(c,n),c);
                %else
                    x_msl(c,n) = NaN;
                %end
            end
        end
        
        % for mhhw
        if z(2501-z_index_mhhw(c,n),c)>= mhhw && z(2501-z_index_mhhw(c,n)-1,c)< mhhw
            
            x_mhhw(c,n) = (mhhw-z(2501-z_index_mhhw(c,n)-1,c))*...
                (x(2501-z_index_mhhw(c,n),c) -  x(2501-z_index_mhhw(c,n)-1,c))/...
                (z(2501-z_index_mhhw(c,n),c) - z(2501-z_index_mhhw(c,n)-1,c))...
                +x(2501-z_index_mhhw(c,n)-1,c);
            
        else if z(2501-z_index_mhhw(c,n),c)<= mhhw && z(2501-z_index_mhhw(c,n)+1,c)> mhhw
                
                x_mhhw(c,n) = (mhhw-z(2501-z_index_mhhw(c,n),c))*...
                    (x(2501-z_index_mhhw(c,n)+1,c) -  x(2501-z_index_mhhw(c,n),c))/...
                    (z(2501-z_index_mhhw(c,n)+1,c) - z(2501-z_index_mhhw(c,n),c))...
                    +x(2501-z_index_mhhw(c,n),c);
                
                % The cross section profile doesn't cross the msl or mhhw
            else %if abs(z_d_mhhw(c,n))< 0.1  % Difference Tolerance < 0.1 meter
                    %x_mhhw(c,n) = x(2501-z_index_mhhw(c,n),c);
                %else
                    x_mhhw(c,n) = NaN;
                %end
            end
        end
        
    end
    
end

%% Visualization Part 1

% Plot graphs on all valid dem dataset (20 in total) at crosssection#n
% for n = 1:25;
n = 10;

% Variable recording data from each dem file for debugging use
x = NaN(2500,17);
z = NaN(2500,17);

figure(1);
% FigH = figure('Position', get(0, 'Screensize'));
grid on; hold on;

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

title(['The Beach Profiles at Transect ',int2str(n)])

set(gca,'FontSize',20)
legend({'Sep 1997','Apr 1998','Sep 2004','Apr 2005','Oct 2005','Mar 2006',...
    'Oct 2006','Apr 2007','Dec 2007','Apr 2008','Oct 2008','Mar 2009',...
    '2009~2011','Sep 2014','Mar 2016','May 2016','Oct 2016'},...
    'Location','SouthEast','Fontsize',10)

xlabel('Cross Shore Distance from 8-meter Isobath (m)')
ylabel('Elevation NAVD88 (m)')

% F    = getframe(FigH);
% imwrite(F.cdata, ['The Beach pofiles at Transect',num2str(n),'a.jpg'], 'jpg')
% close all;
% end

%% Extract single profile graph on dem#m at crosssection#n  Debug Part
% figure(4);grid on; hold on;
% m = 15;   % DEM#m, the x and z data for the single line is x or z(m,1:2500)
% 
% plot(squeeze(packx(m,n,1:packt(m,n))),squeeze(packz(m,n,1:packt(m,n))))
% plot(squeeze(packx(16,n,1:packt(16,n))),squeeze(packz(16,n,1:packt(16,n))))
% 
% xlim([leftb(n),nanmax(x(:,m)+10)])
% 
% title(['Profiles at Cross Section ',int2str(n),' on the DEM dataset ',int2str(m),'&16'])
% xlabel('Cross Shore Distance from -8 Isobath (m)')
% ylabel('Elevation NAVD88 (m)')
% legend('15','16')

%% Visualization Part 2

% % Intervals of transects
% dis_x = fliplr([95,96,106,102,91,110,102,84,103,108,101,98,92,106,92,106,97,...
%     104,111,101,85,114,103,90]);
% x_ab = zeros(1,25);
% 
% % Calculate the along beach distance
% for nl = 2:25
%     x_ab(nl) = x_ab(nl-1) + dis_x(nl-1);
% end
x_ab = 1:25;

%% Calculate the Difference

dmx_msl = NaN(17,25);
dmx_mhhw = NaN(17,25);

% % De-meaned the position by averaged shoreline position at each transect
% mx_msl = nanmean(x_msl);
% mx_mhhw = nanmean(x_mhhw);
% 
% for nl = 1:25
%     dmx_msl(:,nl) = mx_msl(nl) - x_msl(:,nl);
%     dmx_mhhw(:,nl) = mx_mhhw(nl) - x_mhhw(:,nl);
% end


% De-meaned the position by averaged shoreline position at each dataset
mx_msl = nanmean(x_msl,2);
mx_mhhw = nanmean(x_mhhw,2);

for nd = 1:17
    dmx_msl(nd,:) = mx_msl(nd) - x_msl(nd,:);
    dmx_mhhw(nd,:) = mx_mhhw(nd) - x_mhhw(nd,:);
end

% Errorbar
% Preallocate
maxmsl = zeros(1,25); minmsl = zeros(1,25); 
maxmhhw = zeros(1,25); minmhhw = zeros(1,25);

% MSL
for nl = 1:25
maxmsl(nl) = nanmax(dmx_msl(:,nl) - nanmean(dmx_msl(:,nl)));
minmsl(nl) = nanmin(dmx_msl(:,nl) - nanmean(dmx_msl(:,nl)));
end

% MHHW
for nl = 1:25
maxmhhw(nl) = nanmax(dmx_mhhw(:,nl) - nanmean(dmx_mhhw(:,nl)));
minmhhw(nl) = nanmin(dmx_mhhw(:,nl) - nanmean(dmx_mhhw(:,nl)));
end


%% MSL shoreline
figure(2);grid on; hold all;

% De-meaned one
for nd = 1:17
    if nd<=4
    plot(x_ab,dmx_msl(nd,:),'-','LineWidth',1.5);
    else if nd<=8
            plot(x_ab,dmx_msl(nd,:),':','LineWidth',1.5);
        else if nd <=12 
                plot(x_ab,dmx_msl(nd,:),'--','LineWidth',1.5);
            else
                plot(x_ab,dmx_msl(nd,:),'-.','LineWidth',1.5);
            end
        end
    end
end

set(gca,'FontSize',20)
legend({'Sep 1997','Apr 1998','Sep 2004','Apr 2005','Oct 2005','Mar 2006',...
    'Oct 2006','Apr 2007','Dec 2007','Apr 2008','Oct 2008','Mar 2009',...
    '2009~2011','Sep 2014','Mar 2016','May 2016','Oct 2016'},...
    'Location','SouthEast','Fontsize',10)


ylabel('Distance from 8-meter Isobath')
title('Spatially Demeaned Shoreline Position at MSL')
xlabel('Transect Number')
xlim([0,26])
set(gca,'XTick',1:25)
text(-3.5,200,'Accretion','FontSize',15)
text(-3.5,-150,'Erosion','FontSize',15)
% xlabel('Along Beach Distance(m)')
% xlim([-10,2410])


%% MHHW shoreline
figure(3);grid on; hold on;

% De-meaned one
for nd = 1:17
    if nd<=4
    plot(x_ab,dmx_mhhw(nd,:),'-','LineWidth',1.5);
    else if nd<=8
            plot(x_ab,dmx_mhhw(nd,:),':','LineWidth',1.5);
        else if nd <=12 
                plot(x_ab,dmx_mhhw(nd,:),'--','LineWidth',1.5);
            else
                plot(x_ab,dmx_mhhw(nd,:),'-.','LineWidth',1.5);
            end
        end
    end
end

set(gca,'FontSize',20)
legend({'Sep 1997','Apr 1998','Sep 2004','Apr 2005','Oct 2005','Mar 2006',...
    'Oct 2006','Apr 2007','Dec 2007','Apr 2008','Oct 2008','Mar 2009',...
    '2009~2011','Sep 2014','Mar 2016','May 2016','Oct 2016'},...
    'Location','SouthEast','Fontsize',10)


ylabel('Distance from 8-meter Isobath')
title('Spatially Demeaned Shoreline Position at MHHW')
xlabel('Transect Number')
xlim([0,26])
set(gca,'XTick',1:25)
text(-3.5,200,'Accretion','FontSize',15)
text(-3.5,-150,'Erosion','FontSize',15)

% xlabel('Along Beach Distance(m)')
% xlim([-10,2410])


%% Avg MSL and MHHW shoreline
figure(4);grid on; hold all;

% De-meaned Avg
plot(x_ab, nanmean(dmx_msl,1),'-r','LineWidth',1.5);
errorbar(x_ab,nanmean(dmx_msl,1),minmsl,maxmsl,'r')

plot(x_ab, nanmean(dmx_mhhw,1),'-b','LineWidth',1.5);
errorbar(x_ab,nanmean(dmx_mhhw,1),minmhhw,maxmhhw,'b')

set(gca,'FontSize',20)
legend({'Temporally Averaged Shoreline Position at MSL',...
    'Errorbar for Shoreline Position at MSL',...
    'Temporally Averaged Shoreline Position at MHHW',...
    'Errorbar for Shoreline Position at MHHW'},'Fontsize',10,...
    'Location','SouthEast')



title('Spatially Demeaned Shoreline Position')
xlabel('Transect Number')
ylabel('Shoreline Position (meter)')
xlim([0,26])
set(gca,'XTick',1:25)
text(-3.5,250,'Accretion','FontSize',15)
text(-3.3,-200,'Erosion','FontSize',15)
% xlim([-50,2450])
% set (gca,'Xdir','reverse')

