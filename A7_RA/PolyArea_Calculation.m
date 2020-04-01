%% Description
% The script is one of the final scripts CEE 199. 

% A numerical method is used to calculate the volume of the beach and the
% results are compared with that of Surface Volume Tool in Arcgis.

% Linear Regression relationship is used to check the relationships between
% these two method.

clear all; close all; clc

%% Initialization
% Read CSV files
rawdata = csvread('TransectsData_afex_modified.csv',1,1);
Volume_A = csvread('Volume arc.csv',1,1);

% Mean sea level and mean higher high water (Los Angeles Tide Gauge DY1083)
navd = .063;
msl = .861 - navd;
mhhw= 1.674 - navd;

% Preallocate for data, 17 DEMs all have 25 lines, which contain less than
% 1000 data points
packx = NaN(17,25,2500);   % Pack for distance from -8isobath
packz = NaN(17,25,2500);   % Pack for corresponding elevation
packh = zeros(17,25,2500);   % Pack for corresponding height from MSL
packt = zeros(17,25);      % Pack for counting the number of data pairs

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
    
    if rawdata(a,2) < 0 %|| rawdata(a,1) < 700 || rawdata(a,2) > 10
        packz(nd,nl,np) = NaN; % Set NaN
    else
        packz(nd,nl,np) = rawdata(a,2);
    end
    
    % Fix unvalid high elevation point
    if rawdata(a,2) > 10
        packz(nd,nl,np) = (rawdata(a-1, 2) + rawdata(a+1, 2))/2;
        packh(nd,nl,np) = packz(nd,nl,np) - 0.798;
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

%% Calculate the Area Under Curve (Above MSL)
Area = zeros(17,25);

for nd = 1:17
    for nl = 1:25
        for np = 2: packt(nd,nl)
            height_m = (packh(nd,nl,np-1) + packh(nd,nl,np))/2;
            dis_x = packx(nd,nl,np) - packx(nd,nl,np-1);
            Area_sub = dis_x*height_m;
            Area(nd,nl) = Area(nd,nl) + Area_sub;
        end
    end
end

dis_alongb = fliplr([95,96,106,102,91,110,102,84,103,108,101,98,92,106,92,106,97,...
    104,111,101,85,114,103,90]);

Volume_n = zeros(1,17);

for nd = 1:17
    for nl = 2:25
        Area_m = (Area(nd,nl-1) + Area(nd,nl))/2;
        Volume_n(nd) = Volume_n(nd) + Area_m*dis_alongb(nl-1);
    end
end

%% Visualization
y = (Volume_n./2397)'; x = Volume_A./2379;
figure(1);hold on; grid on;
b1 = x\y;
ycalc = b1*x;

plot(x,ycalc,'LineWidth',1.5)
sca = scatter (x,y,50,'filled');
% xlabel('Volume of Method 2: Numerical Method (m^{3}/m)')
% ylabel('Volume of Method 1: Surface Volume Tool (m^{3}/m)')
% plot(ycalc,x,'LineWidth',1.5)
% sca = scatter (y,x,50,'filled');
ylabel('Volume of Method 2: Numerical Method (m^{3}/m)')
xlabel('Volume of Method 1: Surface Volume Tool (m^{3}/m)')


title('Linear Regression Relation Between Two Methods')
legend('Regression Line','Data Pairs','Location','SouthEast')
set(gca,'FontSize',20)
% axis square

Rsq1 = 1 - sum((y - ycalc).^2)/sum((y - mean(y)).^2);
text(380, 395, ['R^{2} = ' num2str(Rsq1)],'FontSize',17);
% text(380, 393, ['y = ',num2str(Rsq1),'x'])