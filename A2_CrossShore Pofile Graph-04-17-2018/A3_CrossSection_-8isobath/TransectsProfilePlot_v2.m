%% Description

%{
Here the accumulative distance value starts from -8m isobath. The profile 
of each transects are plotted but the graphs are not modified, i.e., the
wave affected offshore elevation data is not treated properly.

The shoreline positions at MSL and MHHW are calculated, i.e., find the
distance at which the transect profiles reach MSL or MHHW
%}

clear all; close all; clc

%% Initialization
% Read CSV files
rawdata = csvread('LinePack1_modified.csv',1,1);

% Mean sea level and mean higher high water (Los Angeles Tide Gauge DY1083)
navd = .063;
msl = .861 - navd;
mhhw= 1.674 - navd;

% Preallocate for data, 24 DEMs all have 25 lines, which contain less than
% 2500 data points
packx = NaN(24,25,2500);   % Pack for distance from -8isobath
packz = NaN(24,25,2500);   % Pack for corresponding elevation
packt = zeros(24,25);      % Pack for counting the number of data pairs

% Initialize counting variables for dem, line and data points
nd = 1; nl = 1; np = 1; t = length(rawdata(:,3));

%% Replace invalid value with NaN
% for i = 1: t
%     if rawdata(i,2) == -99
%         rawdata(i,2) = NaN; % Set NaN
%     end
% end
% 
% for i = 2:t-1
%     if ~isnan(rawdata(i,2))
%     if (isnan(rawdata(i-1,2)) && isnan(rawdata(i+1,2)))...
%             || abs(abs(rawdata(i,2))-abs(rawdata(i-1,2)))>5 && abs(abs(rawdata(i,2))-abs(rawdata(i+1,2)))>5
%         rawdata(i,2) = NaN;
%     end
%     end
% end

%% Store data in to data package
for a = 1: t
    
    % Convert data from raw data to datapack
    packx(nd,nl,np) = rawdata(a,1);
    
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
    else
        if rawdata(a+1,3) - rawdata(a,3) ~= 0
            packt(nd,nl) = np-1;
            nd = nd+1;
            nl = 1;
            np = 1;
        end
    end
    end
end


%% Calculate MSL and MHHW distance
% Preallocate for x and z
x = NaN(2500,24);
z = NaN(2500,24);

% Preallocate for recording MSL and MHHW distance
x_msl = NaN(24,25);
x_mhhw = NaN(24,25);

for n = 1:25   % 25 crosssections
    
    for c = 1 : 24   % 24 dem datasets
        x(1:packt(c,n),c) = squeeze(packx(c,n,1:packt(c,n)));
        z(1:packt(c,n),c) = squeeze(packz(c,n,1:packt(c,n)));
        % find index for msl and mhhw
        [z_d_msl, z_index_msl] = min(abs(z - msl));
        [z_d_mhhw, z_index_mhhw] = min(abs(z - mhhw));
        
        % Find corresponding distance
        x_msl(c,n) = x(z_index_msl(c),c) ;
        x_mhhw(c,n) = x(z_index_mhhw(c),c) ;
    end
    
end

x_msl = x_msl';
x_mhhw = x_mhhw';

%% Debug part
% Plot crosssection#n
n = 6;

% Variable recording data from each dem file for debugging use
x = NaN(2500,24);
z = NaN(2500,24);

% plot profiles at different transects
figure(1);grid on; hold on; xlim([700,1200])

for b = 1: 20
x(1:packt(b,n),b) = squeeze(packx(b,n,1:packt(b,n)));
z(1:packt(b,n),b) = squeeze(packz(b,n,1:packt(b,n)));
plot(squeeze(packx(b,n,1:packt(b,n))),squeeze(packz(b,n,1:packt(b,n))))
end

legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'...
    ,'16','17','18','19','20')

[z_d_msl, z_index_msl] = min(abs(z - msl));
[z_d_mhhw, z_index_mhhw] = min(abs(z - mhhw));

% plot
figure(2);grid on; hold on; xlim([700,1200])
sb = 20;
plot(squeeze(packx(sb,n,1:packt(sb,n))),squeeze(packz(sb,n,1:packt(sb,n))))
