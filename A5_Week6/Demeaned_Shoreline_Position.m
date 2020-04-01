%% Description
% The code is for a CEE 199 assignment

% 24 DEM are used in arcgis to extrac cross shore profile of the Sunset
% Beach. The last 4 DEM data are considered as invalid and thus only 20 DEM
% are used to plot

%This Script plot the mean shoreline position at different date

% The shoreline positions are demeaned

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

%% Calculate MSL and MHHW distance

% Preallocate for x and z (Set an oversized matrix, because different dem
% have different data set size.

x = NaN(2500,24);
z = NaN(2500,24);

% Preallocate for recording MSL and MHHW distance
x_msl = NaN(24,25);
x_mhhw = NaN(24,25);
z_d_msl = NaN(24,25); z_index_msl = NaN(24,25);
z_d_mhhw = NaN(24,25); z_index_mhhw = NaN(24,25);

for n = 1:25   % 25 crosssections
    
    for c = 1 : 20   % 24 dem datasets
        x(1:packt(c,n),c) = squeeze(packx(c,n,1:packt(c,n)));
        z(1:packt(c,n),c) = squeeze(packz(c,n,1:packt(c,n)));
        
        
        % find index for msl and mhhw
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
x_dif_msl = NaN(15,25);
x_dif_mhhw = NaN(15,25);
m_msl =  nanmean(x_msl(1:15,:));
m_mhhw = nanmean(x_mhhw(1:15,:));

for nl = 1:25
    x_dif_msl(:,nl) = m_msl(nl) - x_msl(1:15,nl);
    x_dif_mhhw(:,nl) = m_mhhw(nl) - x_mhhw(1:15,nl);
end

m_dif_msl = mean(x_dif_msl,2);
m_dif_mhhw = mean(x_dif_mhhw,2);

% Dataset 13 "2009-2011 CA Coastal Conservancy Coastal Lidar Project_DEM"
% don't have specific collecting date
m_dif_msl(13) = NaN;
m_dif_mhhw(13) = NaN;

%% Visualization
Datestring = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
    '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';'12/2009';...
    '09/2014';'05/2016'];
% Datestring = ['04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
%     '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';...
%     '05/2016'];
formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);
ref = zeros(15,1);

% For MSL
figure(1); hold all;grid on;
plot(ta,ref)
plot(ta,m_dif_msl)
% plot(ta,[m_dif_msl(2:12);m_dif_msl(15)])
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Mean Shore Line Position (m)')
title ('MSL')

% For MHHW
figure(2); hold all;grid on;
plot(ta,ref)
plot(ta,m_dif_mhhw)
% plot(ta,[m_dif_msl(2:12);m_dif_msl(15)])
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Mean Shore Line Position (m)')
title ('MHHW')
%% Visualization partion 1
% 
% % Plot graphs on all valid dem dataset (20 in total) at crosssection#n
% % for n = 1:25;
% n = 6;
% 
% % Variable recording data from each dem file for debugging use
% x = NaN(2500,24);
% z = NaN(2500,24);
% 
% figure(1);grid on; hold on;
% 
% for b = 1: 15
%     x(1:packt(b,n),b) = squeeze(packx(b,n,1:packt(b,n)));
%     z(1:packt(b,n),b) = squeeze(packz(b,n,1:packt(b,n)));
%     if b<= 5
%         plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
%             ,'-','LineWidth',2)
%     else if b<=10
%             plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
%                 ,':','LineWidth',2)
%         else if b<=15
%                 plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
%                     ,'-.','LineWidth',2)
%             else
%                 plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n)))...
%                     ,'--','LineWidth',2)
%             end
%         end
%     end
% end
% 
% % Extract the maximum distance from the farthest point of the longest
% % profile to the -8 isobath
% 
% leftb = [720,750,800,850,900,900,850,850,820,850,850,850,750,920,920,920,930,980,980,950,920,900,900,880,850];
% xlim([leftb(n),nanmax(nanmax(squeeze(packx(:,n,:))))+10])
% 
% title(['Profiles at Cross Section ',int2str(n)])
% legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'...
%     ,'Location','SouthEast')
% xlabel('Distance from -8 Isobath (m)')
% ylabel('Elevation (m)')
% 
% % % % Debug part
% % pause (10)
% % close all
% % end
% 
%% Extract single profile graph on dem#m at crosssection#n
% figure(2);grid on; hold on;
% m = 14;   % DEM#m, the x and z data for the single line is x or z(m,1:2500)
% 
% plot(squeeze(packx(m,n,1:packt(m,n))),squeeze(packz(m,n,1:packt(m,n))))
% xlim([leftb(n),nanmax(x(:,b)+10)])
% 
% title(['Profiles at Cross Section ',int2str(n),' on the DEM dataset ',int2str(m)])
% xlabel('Distance from -8 Isobath (m)')
% ylabel('Elevation (m)')

