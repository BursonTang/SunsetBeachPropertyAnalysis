%% Description
% The is the tide series generated with observed water levels at La Jolla,
% Los Angeles and San Diego from 12/01/2015 00:00 to 2015/12/31 23:59 GTM.
% NOAA : https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels
clear all;close all;clc;

%% Read file and Extract data
% Read CSV files
f1 = csvread('CO-OPS__LaJolla__hr.csv',1,1);
f2 = csvread('CO-OPS__LosAngeles__hr.csv',1,1);
f3 = csvread('CO-OPS__SanDiego__hr.csv',1,1);

% Extract wave data
wave_lj = f1(:,1).';
wave_la = f2(:,1).';
wave_sd = f3(:,1).';

%% Plot 
% generate datetime series
t = datetime(2015,12,1) + hours(0:743);


figure(1);hold on;grid on;
plot(t,wave_lj,'r','LineWidth',1.5)
plot(t,wave_la,'g','LineWidth',1.5)
plot(t,wave_sd,'b','LineWidth',1.5)

% Set labels
xlabel('Date (UTC)')
ylabel('Height above MLLW(m)')
title('Verfied Water Level at Three Stations (From 2015/12/01 00:00GTM to 2015/12/31 23:59GTM)')
legend('La Jolla Station','Los Angeles Station','San Diego Station','Location','SouthWest')
gh = gca;
gh.FontSize = 20;