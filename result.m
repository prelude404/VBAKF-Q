function [error_xyz,error] = result(x_esti, gtd, imu, uwb, t, fig_num)
% clf;
figure(fig_num)
set(gcf,'Position',[100,20,600,600]);
subplot(4,1,1)
plot3(x_esti(1,:),x_esti(2,:),x_esti(3,:),'b--',gtd(1,:),gtd(2,:),gtd(3,:),'b-','linewidth',1);
legend('x_{esti}','x_{gtd}');
ylabel('Trajectory');
title('Result');
subplot(4,1,2)
plot(t,gtd(1,:),'r-',t,gtd(2,:),'m-',t,gtd(3,:),'b-','linewidth',1);
hold on
plot(t,x_esti(1,:),'r--',t,x_esti(2,:),'m--',t,x_esti(3,:),'b--','linewidth',1);
h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','Location','northwest');
ylabel('Position');
set(h1,'Orientation','horizon','Box','on');
subplot(4,1,3)
plot(t,gtd(4,:),'r-',t,gtd(5,:),'m-',t,gtd(6,:),'b-','linewidth',1);
hold on
plot(t,x_esti(4,:),'r--',t,x_esti(5,:),'m--',t,x_esti(6,:),'b--','linewidth',1);
h2 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','Location','northwest');
ylabel('Velocity')
set(h2,'Orientation','horizon','Box','on');
subplot(4,1,4)
error = gtd - x_esti;
error_norm = sqrt(error(1,:).^2 + error(2,:).^2 + error(3,:).^2);
plot(t,error_norm,'k','linewidth',1);
xlabel('Time')
ylabel('Estimation Error')
legend('error')
error_xyz = zeros(1,4);
error_xyz(1) = sqrt(mean(error(1,:).^2));
error_xyz(2) = sqrt(mean(error(2,:).^2));
error_xyz(3) = sqrt(mean(error(3,:).^2));
error_xyz(4) = sqrt(mean(error(1,:).^2 + error(2,:).^2 + error(3,:).^2));