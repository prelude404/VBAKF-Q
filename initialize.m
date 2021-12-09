function [imu_noise, uwb_noise, K, dt, t] = initialize()
% clc;
clear;
dt = 1/25;
K = 300/dt;
t = 0:dt:K*dt;
t(K+1) = [];
T = 2;

sigma_u = diag([0.01,0.01,0.01]);
sigma_y = 0.01;
wave_imu = [1;1;1]*(0.55 + 0.45 * sin(2*pi/T*t));
wave_uwb = 1 + 0.05 * sin(2*pi/(2*T)*t);

imu_noise = wave_imu .* (sqrt(sigma_u)*randn(3,K));
uwb_noise = wave_uwb .* (sqrt(sigma_y)*randn(1,K));

