function [gtd, u, y, imu, uwb] = curve(imu_noise, uwb_noise, t)
x0 = [50;50;50];
T = 20;

p_curve = [x0(1)+200*(cos(2*pi/T*t)-1);
           x0(2)+200*sin(1.5*pi/T*t);
           x0(3)+200*sin(1*pi/T*t)];
v_curve = [-200*(2*pi/T)*sin(2*pi/T*t);
           200*(1.5*pi/T)*cos(1.5*pi/T*t);
           200*(1*pi/T)*cos(1*pi/T*t)];
u_curve = [-200*((2*pi/T)^2)*cos(2*pi/T*t);
           -200*((1.5*pi/T)^2)*sin(1.5*pi/T*t);
           -200*((1*pi/T)^2)*sin(1*pi/T*t)];

gtd(1:3,:) = p_curve;
gtd(4:6,:) = v_curve;
gtd(7,:) = p_curve(:,1)' * v_curve(:,1);
gtd(8,:) = v_curve(:,1)' * v_curve(:,1);

u = u_curve;
y = sqrt(p_curve(1,:).^2 + p_curve(2,:).^2 + p_curve(3,:).^2);
imu = u + imu_noise;
uwb = y + uwb_noise;

