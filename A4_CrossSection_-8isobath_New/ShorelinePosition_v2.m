%% Description
%{
Here the accumulative distance value starts from -8m isobath. The profile 
of each transects are plotted but the graphs are not modified, i.e., the
wave affected offshore elevation data is not treated properly.

The shoreline positions at MSL and MHHW are calculated, i.e., find the
distance at which the transect profiles reach MSL or MHHW

At some of the date, the dem collected doesn't reach MSL or MHHW and thus
eleminated in the plot, so only 16 dem data are used in this Script
%}

clear all; close all; clc

%% For MSL
%Initialization
% Read CSV files
rawdata = csvread('MSL difference.csv',1,1);

% Delete invalid data (data with confusing collecting data)
rawdata_d = [rawdata(1:11, :); rawdata(16,:)];

Datestring = ['04/1998';'09/2004';'04/2005';'10/2005';'03/2006';'10/2006';...
    '04/2007';'12/2007';'04/2008';'10/2008';'03/2009';'05/2016'];
formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);

%% Plot Individually
% For Cross Section - n
figure(1); hold on;grid on;
n = 1;
ref = zeros(12,1);
plot(ta,ref)
plot(ta,rawdata_d(:,n))
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title (['MSL Position at transect',num2str(n)])

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
rawdata = csvread('MHHW difference.csv',1,1);

% Delete invalid data (data with confusing collecting data)
rawdata_d = [rawdata(1:11, :); rawdata(16,:)];

Datestring = ['04/1998';'09/2004';'04/2005';'10/2005';'03/2006';'10/2006';...
    '04/2007';'12/2007';'04/2008';'10/2008';'03/2009';'05/2016'];
formatin = 'mm/yyyy';
ta = datenum(Datestring,formatin);

%% Plot Individually
% For Cross Section - n
figure(2); hold on;grid on;
n = 1;
ref = zeros(12,1);
plot(ta,ref)
plot(ta,rawdata_d(:,n))
dateFormat = 'mm/yyyy';
datetick('x',dateFormat)
xlabel('Date')
ylabel('Distance from 8m Isobath (m)')
title (['MHHW Position at transect',num2str(n)])

%% Subplot
% for i = 1:25
%     subplot(5,5,i);hold on;
%     plot(ta,rawdata_d(:,i),'-')
%     dateFormat = 'mm/yyyy';
%     datetick('x',dateFormat)
%     grid on; 
% end
