function [x_kf] = kf(gtd, imu, uwb, t)
% kf: The classical Kalman filter
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r

%% Data Preparation: y,u,K,dt
x0 = gtd(1:6,1);
K = length(t);
dt = t(2) - t(1);
u = imu;
v = zeros(3,K);
d = zeros(3,K);
y = zeros(1,K);
y0 = 0.5 * (uwb(1).^2);
y1 = 0.5 * (uwb.*uwb);
for i = 1:K
    v(1,i) = dt * trapz(u(1,1:i));
    v(2,i) = dt * trapz(u(2,1:i));
    v(3,i) = dt * trapz(u(3,1:i));
    d(1,i) = dt * trapz(v(1,1:i));
    d(2,i) = dt * trapz(v(2,1:i));
    d(3,i) = dt * trapz(v(3,1:i));

    y(i) = y1(i) - y0 + 0.5 * (d(1,i)^2 + d(2,i)^2 + d(3,i)^2);
end

%% x_k = A x_k-1 + B u_k
A = [1,0,0,dt,0,0,0,0;
     0,1,0,0,dt,0,0,0;
     0,0,1,0,0,dt,0,0;
     0,0,0, 1,0,0,0,0;
     0,0,0, 0,1,0,0,0;
     0,0,0, 0,0,1,0,0;
     0,0,0, 0,0,0,1,0;
     0,0,0, 0,0,0,0,1];

B = [0.5*dt^2,0,0;
     0,0.5*dt^2,0;
     0,0,0.5*dt^2;
     dt,   0,   0;
     0,   dt,   0;
     0,    0,  dt;
     0,    0,   0;
     0,    0,   0];

%% Initialization
x0 = x0 + [0;0;0;0;0;0]; % add initial error
x_esti = zeros(8,K);
x_esti(:,1) = [x0(1:3,1); x0(4:6,1); x0(1:3,1)'*x0(4:6,1); x0(4:6,1)'*x0(4:6,1)];
x_pre = zeros(8,K);
x_pre(:,1) = x_esti(:,1);

%% Parameter
P = diag([1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-2,1e-1]);
Q = diag([1e-1,1e-1,1e-1,1e-7,1e-7,1e-7,1e-2,1e-4]);
R = 1e-2;

%% KF
for i = 2:K
    x_pre(:,i) = A * x_esti(:,i-1) + B * u(:,i);
    P = A * P * A'+Q;
    H = [d(1,i),d(2,i),d(3,i),0,0,0,((i-1)*dt),0.5*(((i-1)*dt)^2)];
    K = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = (eye(8)-K*H)*P;
end
x_kf = x_esti;


