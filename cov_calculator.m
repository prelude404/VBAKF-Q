function [Q,R] = cov_calculator(imu, uwb, gtd, t)
% gtd_y
K = length(t);
dt = t(2) - t(1);
gtd_v = zeros(3,K);
gtd_d = zeros(3,K);
gtd_y = zeros(1,K);
gtd_y0 = 0.5 * sum(gtd(1:3,1).^2);
gtd_y1 = 0.5 * sum(gtd(1:3,:).^2);
for i = 1:K
    gtd_v(1:3,i) = gtd(4:6,i) - gtd(4:6,1);
    gtd_d(1,i) = dt * trapz(gtd_v(1,1:i));
    gtd_d(2,i) = dt * trapz(gtd_v(2,1:i));
    gtd_d(3,i) = dt * trapz(gtd_v(3,1:i));

    gtd_y(i) = gtd_y1(i) - gtd_y0 + 0.5 * (gtd_d(1,i)^2 + gtd_d(2,i)^2 + gtd_d(3,i)^2);
end

% y
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

% gtd_u
gtd_u = zeros(3,K);
for i = 2:K
    gtd_u(1:3,i) = (gtd(4:6,i) - gtd(4:6,i))./dt;
end

% calculate Q
B = [0.5*dt^2,0,0;
     0,0.5*dt^2,0;
     0,0,0.5*dt^2;
     dt,   0,   0;
     0,   dt,   0;
     0,    0,  dt;
     0,    0,   0;
     0,    0,   0];
error_u = imu - gtd_u;
Q0 = mean(error_u'.^2);
Q = B * diag(Q0) * B';

% calculate R
error_y = y - gtd_y;
R = mean(error_y.^2);

