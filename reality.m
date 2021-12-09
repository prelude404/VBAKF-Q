function [gtd, imu, uwb, K, dt,t]=reality()
% This function deal with data bag in real flight

%% load data
data = load('data/bag1.mat');
% data = load('data/bag2.mat');
% data = load('data/bag3.mat');
% data = load('data/bag4.mat');
% data = load('data/bag5.mat');

starti = 175;endi = 605;   %bag1
% starti = 175;endi = 735;   %bag2
% starti = 115;endi = 740;   %bag3
% starti = 210;endi = 1140;  %bag4
% starti =  70;endi = 1030;  %bag5

%% store data
% time = data.time;       % time
% gtd = data.gtd';        % position from vicon
% gtd(4:6,:) = data.vel'; % velocity from vicon
% uwb = data.uwb;         % uwb distance
% imu = data.imu';        % imu data

t = data.time(starti:endi);
t = t - t(1);
dt = t(2) - t(1);
K = length(t);
uwb = data.uwb(:,starti:endi);
imu = data.imu(starti:endi,:)';
% imu(3,:) = imu(3,:) - 9.7964;
imu(1,:) = imu(1,:) - 0.11; % eliminate the bias of axis X
imu(2,:) = imu(2,:) + 0.01; % eliminate the bias of axis Y
imu(3,:) = imu(3,:) - 10.08; % eliminate the bias of axis Z

gtd(1:3,:) = data.gtd(starti:endi,:)';
gtd(4:6,:) = data.vel(starti:endi,:)';
gtd(7,:) = gtd(1:3,1)' * gtd(4:6,1);
gtd(8,:) = gtd(4:6,1)' * gtd(4:6,1);



