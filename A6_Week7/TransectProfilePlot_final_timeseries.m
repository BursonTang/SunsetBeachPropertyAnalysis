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

%% Calculate the Difference
dmt_msl = NaN(17,25);
dmt_mhhw = NaN(17,25);
m_msl =  nanmean(x_msl);
m_mhhw = nanmean(x_mhhw);

% for nd = 1:17
%     x_dif_msl(nd,:) = m_msl(nd) - x_msl(nd,:);
%     x_dif_mhhw(nd,:) = m_mhhw(nd) - x_mhhw(nd,:);
% end

for nl = 1:25
    dmt_msl(:,nl) = m_msl(nl) - x_msl(:,nl);
    dmt_mhhw(:,nl) = m_mhhw(nl) - x_mhhw(:,nl);
end

m_dif_msl = nanmean(dmt_msl,2);
m_dif_mhhw = nanmean(dmt_mhhw,2);

% Dataset 13 "2009-2011 CA Coastal Conservancy Coastal Lidar Project_DEM"
% don't have specific collecting date
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
maxmsl = zeros(17,1); minmsl = zeros(17,1); 
maxmhhw = zeros(17,1); minmhhw = zeros(17,1);

% MSL
for nd = 1:17
maxmsl(nd) = nanmax( dmt_msl(nd,:) - nanmean(dmt_msl(nd,:)));
minmsl(nd) = nanmin( dmt_msl(nd,:) - nanmean(dmt_msl(nd,:)));
end

% MHHW
for nd = 1:17
maxmhhw(nd) = nanmax( dmt_mhhw(nd,:) - nanmean(dmt_mhhw(nd,:)));
minmhhw(nd) = nanmin( dmt_mhhw(nd,:) - nanmean(dmt_mhhw(nd,:)));
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
    '09/2014';'03/2016';'05/2016';'10/2016'];

formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);
ta(1) = 729750;

% In order to draw a reference line
ref = zeros(16,1);

% For MSL
figure(1); hold all;grid on;
%plot(ta,ref)
plot(ta,[m_dif_msl(1:12);m_dif_msl(14:17)],'-r','LineWidth',1.5)
errorbar(ta,[m_dif_msl(1:12);m_dif_msl(14:17)],[minmsl(1:12);minmsl(14:17)],...
    [maxmsl(1:12);maxmsl(14:17)],'r')
% plot(ta,[m_dif_msl(2:12);m_dif_msl(15)])
% dateFormat = 'mm/yyyy';
% datetick('x',dateFormat)
% xlabel('Date')
% ylabel('Mean Shore Line Position (m)')
% title ('MSL Position Timeseries')

% For MHHW
%figure(2); hold all;grid on;
%plot(ta,ref)
plot(ta,[m_dif_mhhw(1:12);m_dif_mhhw(14:17)],'-b','LineWidth',1.5)
errorbar(ta,[m_dif_mhhw(1:12);m_dif_mhhw(14:17)],[minmhhw(1:12);minmhhw(14:17)],...
    [maxmhhw(1:12);maxmhhw(14:17)],'b')
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Shoreline Position (meter)')
set(gca,'FontSize',20)
title ('Temporally Demeaned Shoreline Position')
text(728700,80,'Accretion','FontSize',15)
text(728750,-60,'Erosion','FontSize',15)
legend({'Spatially Averaged Shoreline Position at MSL',...
    'Errorbar for Shoreline Position at MSL',...
    'Spatially Averaged Shoreline Position at MHHW',...
    'Errorbar for Shoreline Position at MHHW'},'Location','SouthWest','Fontsize',10)
