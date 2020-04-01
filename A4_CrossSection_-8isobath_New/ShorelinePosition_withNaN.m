%% Description
%{
Here the accumulative distance value starts from -8m isobath. The profile 
of each transects are plotted but the graphs are not modified, i.e., the
wave affected offshore elevation data is not treated properly.

The shoreline positions at MSL and MHHW are calculated, i.e., find the
distance at which the transect profiles reach MSL or MHHW

Although at some of the date, the dem collected doesn't reach MSL or MHHW,
they are used in this Script, totally 20 dem, and thus 20 date
%}


clear all; close all; clc

%% For MSL
%Initialization
% Read CSV files
rawdata = csvread('MSL difference_NaN.csv',1,1);

% Delete invalid data (data with confusing collecting data)
rawdata_d = [rawdata(1:12, :); rawdata(14, :); rawdata(19:20,:)];

Datestring = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
    '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';'12/2009';...
    '09/2014';'05/2016'];
formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);

%% Plot Individually

% For Cross Section - n
figure(1); hold on;grid on;
n = 10;
ref = zeros(15,1);
plot(ta,ref) % Plot a reference line
plot(ta,rawdata_d(:,n))

dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title (['MSL Position at transect',num2str(n)])

% Mean location
figure(2); hold on; grid on;
msl_avg_d = nansum(rawdata_d,2)/25;

plot(ta,ref)
plot(ta,msl_avg_d)
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title ('Mean MSL Position')
%% Subplot
% for i = 1:25
%     subplot(5,5,i);hold on;
%     plot(ta,rawdata_d(:,i),'-')
%     dateFormat = 'mm/yyyy';
%     datetick('x',dateFormat)
%     grid on; 
% end

%% For MHHW
%Initialization
% Read CSV files
rawdata = csvread('MHHW difference_NaN.csv',1,1);

% Delete invalid data (data with confusing collecting data)
rawdata_d = [rawdata(1:12, :); rawdata(14, :); rawdata(19:20,:)];

Datestring = ['09/1997';'04/1998';'09/2004';'04/2005';'10/2005';'03/2006';...
    '10/2006';'04/2007';'12/2007';'04/2008';'10/2008';'03/2009';'12/2009';...
    '09/2014';'05/2016'];
formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);

%% Plot Individually
% For Cross Section - n
figure(3); hold on; grid on;
% n = 10;
ref = zeros(15,1);
plot(ta,ref)
plot(ta,rawdata_d(:,n))
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title (['MHHW Position at transect',num2str(n)])

% Mean location
figure(4); hold on; grid on;
mhhw_avg_d = nansum(rawdata_d')/25;

plot(ta,ref)
plot(ta,mhhw_avg_d)
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title ('Mean MHHW Position')
%% Subplot
% for i = 1:25
%     subplot(5,5,i);hold on;
%     plot(ta,rawdata_d(:,i),'-')
%     dateFormat = 'mm/yyyy';
%     datetick('x',dateFormat)
%     grid on; 
% end
