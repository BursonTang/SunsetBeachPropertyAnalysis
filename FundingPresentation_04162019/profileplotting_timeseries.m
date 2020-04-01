%% Description
% The script is one of the final scripts CEE 199. 

% This script plots the demeaned shoreline position at different date. The
% shoreline positions are are demeaned by averaged shoreline 
% position at each transect

clear all; close all; clc

%% Initialization
% Read CSV files
rawdata = csvread('TransectsData_afex_modified.csv',1,1);

% Mean sea level and mean higher high water (Los Angeles Tide Gauge DY1083)
navd = .063;
msl = .861 - navd;
mhhw= 1.674 - navd;

% Preallocate for data, 24 DEMs all have 25 lines, which contain less than
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
for c = 1:17   % 24 dem datasets
    
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

atv_msl = [1040.8016, 1035.9031, 1032.1551, 1025.4141, 1034.9884, 1050.38, ...
    1063.242, 1078.6053, 1041.6894, 1012.7816, 1019.0022, 989.9192,...
    856.7471, 904.4320, 919.9711, 931.2243, 898.8577, 916.2077,...
    931.4211, 972.1354, 966.2935, 915.1720, 857.4641, 810.3354, 782.5382];

atv_mhhw = [1048.8865, 1042.7998, 1036.49, 1030.9627, 1040.2638, 1055.222, 1069.404,...
    1084.896, 1047.4506, 1018.6423, 1024.3552, 996.7297, 863.5496,...
    910.8627, 925.903, 937.6273, 905.5505, 923.3011, 938.1243,...
    978.9436, 975.629, 922.3396, 865.3639, 817.6038, 792.7410];
%% Calculate the Difference
% add new dataset
x_msl = [x_msl;atv_msl];
x_mhhw = [x_mhhw;atv_mhhw];

% preallocate for difference of shoreline position
dmt_msl = NaN(18,25); dmt_mhhw = NaN(18,25);

% calculate temporally averaged shoreline position
m_msl =  nanmean(x_msl); m_mhhw = nanmean(x_mhhw);

% for nd = 1:17
%     x_dif_msl(nd,:) = m_msl(nd) - x_msl(nd,:);
%     x_dif_mhhw(nd,:) = m_mhhw(nd) - x_mhhw(nd,:);
% end

% Calculate the erosion and accretion of the beach (time averaged SLP - individual SLP)
for nl = 1:25
    dmt_msl(:,nl) = m_msl(nl) - x_msl(:,nl);
    dmt_mhhw(:,nl) = m_mhhw(nl) - x_mhhw(:,nl);
end

% calculate the mean erosion and accretion at each time period, i.e.,
% average erosion and accretion along the beach
m_dif_msl = nanmean(dmt_msl,2);
m_dif_mhhw = nanmean(dmt_mhhw,2);

% Dataset 13 "2009-2011 CA Coastal Conservancy Coastal Lidar Project_DEM"
% don't have specific collecting date, and thus cannot be plotted on
% timeseeries, set as NaN
m_dif_msl(13) = NaN;
m_dif_mhhw(13) = NaN;

% Ignore datasets with too many NaN values
% for i = 1:17
%     
%     if length( find(isnan(x_msl(i,:))) ) >10
%         m_dif_msl(i) = NaN;
%     end
%     
%     if length( find(isnan(x_mhhw(i,:))) ) >10
%         m_dif_mhhw(i) = NaN;
%     end
%     
% end


% Errorbar
% Preallocate
Umsl = zeros(18,1); Lmsl = zeros(18,1); 
Umhhw = zeros(18,1); Lmhhw = zeros(18,1);

% MSL
for nd = 1:18
Umsl(nd) = nanmax( dmt_msl(nd,:) - nanmean(dmt_msl(nd,:)));
Lmsl(nd) = nanmin( dmt_msl(nd,:) - nanmean(dmt_msl(nd,:)));
end

% MHHW
for nd = 1:18
Umhhw(nd) = nanmax( dmt_mhhw(nd,:) - nanmean(dmt_mhhw(nd,:)));
Lmhhw(nd) = nanmin( dmt_mhhw(nd,:) - nanmean(dmt_mhhw(nd,:)));
end

% % Errorbar
% % MSL
% maxmsl = nanmax(dmt_msl,[],2);
% minmsl = nanmin(dmt_msl,[],2);
% 
% % MHHW
% maxmhhw = nanmax(dmt_mhhw,[],2);
% minmhhw = nanmin(dmt_mhhw,[],2);
%% Visualization
% The collceted date of new 17 datasets(week7),
Datestring = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
    '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';...
    '09/2014';'03/2016';'05/2016';'10/2016';'11/2017'];


formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);
Datestring2 = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
    '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';...
    '09/2014';'03/2016';'05/2016';'11/2017'];%;'05/2016';
ta2 = datenum(Datestring2,formatin);
ta(1) = 729750;

% In order to draw a reference line
ref = zeros(16,1);

% Plot for MSL
figure(1); hold all;grid on;
%plot(ta,ref)
plot(ta2,[m_dif_msl(1:12);m_dif_msl(14:16);m_dif_msl(18)],'-b','LineWidth',1.5)
errorbar(ta,[m_dif_msl(1:12);m_dif_msl(14:18)],[Lmsl(1:12);Lmsl(14:18)],...
    [Umsl(1:12);Umsl(14:18)],'b')

% Plot for MHHW
plot(ta,[m_dif_mhhw(1:12);m_dif_mhhw(14:18)],'-r','LineWidth',1.5)
errorbar(ta,[m_dif_mhhw(1:12);m_dif_mhhw(14:18)],[Lmhhw(1:12);Lmhhw(14:18)],...
    [Umhhw(1:12);Umhhw(14:18)],'r')

%Plot points for beach nourishment projects
yposition(1) = m_dif_msl(2)+(m_dif_msl(3)-m_dif_msl(2))/(ta(3)-ta(2)) * (datenum('2002','yyyy')-ta(2));
yposition(2) = m_dif_msl(12);
% plot([datenum('2002','yyyy'),datenum('2009','yyyy')],)
scatter([datenum('2002','yyyy'),ta(12)],yposition(:),'pk','LineWidth',5)

% Graph Setting
% Set the Xaxis as date from 2002 to 2017
gcahand = gca;
gcahand.XTick = datenum( datetime(2002,1,2) + years(0:17) );
gcahand.XLim = [datenum('11/2001','mm/yyyy'),datenum('01/02/2018','mm/dd/yyyy')];
dateFormat = 'yyyy';%dateFormat = 'mm/yyyy';
datetick('x',dateFormat,'keeplimits','keepticks')
gcahand.XTickLabelRotation = 45; % Rotate the Xtick label

% set titles
xlabel('Year')
ylabel('Shoreline Position change (m)')
title ('Average Shoreline Position Change')
% set fontsize
gcahand.FontSize = 20;

% Legend and explanation text
text(datenum('11/2001','mm/yyyy')-620,60,'Accretion','FontSize',15)
text(datenum('11/2001','mm/yyyy')-600,-60,'Erosion','FontSize',15)
legend({'Spatially Averaged Shoreline Position at MSL',...
    'Min/Max Average MSL Shoreline Position',...
    'Spatially Averaged Shoreline Position at MHHW',...
    'Min/Max Average MHHW Shoreline Position'...
    'Beach Nourishment Project'},'Location','SouthWest','Fontsize',10)

%% plot only MSL
% figure(2); hold all;grid on;
% %plot(ta,ref)
% plot(ta,[m_dif_msl(1:12);m_dif_msl(14:17)],'-b','LineWidth',1.5)
% errorbar(ta,[m_dif_msl(1:12);m_dif_msl(14:17)],[Lmsl(1:12);Lmsl(14:17)],...
%     [Umsl(1:12);Umsl(14:17)],'b')
% 
% %Plot points for beach nourishment projects
% yposition(1) = m_dif_msl(2)+(m_dif_msl(3)-m_dif_msl(2))/(ta(3)-ta(2)) * (datenum('2002','yyyy')-ta(2));
% yposition(2) = m_dif_msl(12);
% % plot([datenum('2002','yyyy'),datenum('2009','yyyy')],)
% scatter([datenum('2002','yyyy'),ta(12)],yposition(:),'pr','LineWidth',5)
% 
% % Graph Setting
% % Set the Xaxis as date from 2002 to 2017
% gcahand = gca;
% gcahand.XTick = datenum( datetime(2002,1,2) + years(0:16) );
% gcahand.XLim = [datenum('11/2001','mm/yyyy'),datenum('01/02/2017','mm/dd/yyyy')];
% dateFormat = 'yyyy';%dateFormat = 'mm/yyyy';
% datetick('x',dateFormat,'keeplimits','keepticks')
% gcahand.XTickLabelRotation = 45; % Rotate the Xtick label
% 
% % set titles
% xlabel('Year')
% ylabel('Shoreline Position change (m)')
% title ('Average Shoreline Position Change')
% % set fontsize
% gcahand.FontSize = 20;
% 
% % Legend and explanation text
% text(datenum('11/2001','mm/yyyy')-620,60,'Accretion','FontSize',15)
% text(datenum('11/2001','mm/yyyy')-600,-60,'Erosion','FontSize',15)
% legend({'Spatially Averaged Shoreline Position at MSL',...
%     'Deviation of Shoreline Position at MSL',...
%     'Beach Nourishment Project'},'Location','SouthWest','Fontsize',10)