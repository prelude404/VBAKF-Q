close all;
[imu_noise, uwb_noise, K, dt, t] = initialize();
[gtd, u, y, imu, uwb] = curve(imu_noise, uwb_noise, t);
[x_esti, x_predict] = vbakf_q(gtd, imu, uwb, t);
[error_xyz, error] = result(x_esti, gtd, imu, uwb, t, 1);
disp(['VBAKF-Q ErrorX: ',num2str(error_xyz(1)),'  ErrorY: ',num2str(error_xyz(2)),'  ErrorZ: ',num2str(error_xyz(3)),'  ErrorTotal: ',num2str(error_xyz(4))]);

% [x_esti1, x_predict1,d] = akf_mhe(gtd, imu, uwb, t);
% [error_xyz1, error1] = result(x_esti1, gtd, imu, uwb, t, 2);
% disp(['AKF-MHE ErrorX: ',num2str(error_xyz1(1)),'  ErrorY: ',num2str(error_xyz1(2)),'  ErrorZ: ',num2str(error_xyz1(3)),'  ErrorTotal: ',num2str(error_xyz1(4))]);

%% compare with classical KF
% [x_kf] = kf(gtd, imu, uwb, t);
% [error0_xyz,error0] = result(x_kf, gtd, imu, uwb, t, 2);
% disp(['ErrorX0: ',num2str(error0_xyz(1)),'  ErrorY0: ',num2str(error0_xyz(2)),'  ErrorZ0: ',num2str(error0_xyz(3)),'  ErrorTotal0: ',num2str(error0_xyz(4))]);

%% deal with real data
% [gtd, imu, uwb, K, dt,t] = reality();
% [x_esti, x_predict] = vbakf_q(gtd, imu, uwb, t);
% 
% [error_xyz,error] = result(x_esti, gtd, imu, uwb, t, 1);
% disp(['ErrorX: ',num2str(error_xyz(1)),'  ErrorY: ',num2str(error_xyz(2)),'  ErrorZ: ',num2str(error_xyz(3)),'  ErrorTotal: ',num2str(error_xyz(4))]);

%% Output
% figure(3)
% plot3(gtd(1,:),gtd(2,:),gtd(3,:),'b-','linewidth',1);
% xlabel('x','FontName','Times New Roman','FontSize',16);
% ylabel('y','FontName','Times New Roman','FontSize',16);
% zlabel('z','FontName','Times New Roman','FontSize',16);
% title('Trajectory','FontName','Times New Roman','FontSize',16);
% grid on;
% 
% figure(4)
% axis = 3;
% plot(t,gtd(axis,:),'r-',t,x_esti(axis,:),'m-',t,x_predict(axis,:),'b-','linewidth',1);
% hold on
% h1 = legend('gtd','esti','predict','FontName','Times New Roman','FontSize',12);
% xlabel('Time','FontName','Times New Roman','FontSize',16);
% ylabel('Position','FontName','Times New Roman','FontSize',16);
% set(h1,'Orientation','horizon','Box','on');
% title('Estimation Result','FontName','Times New Roman','FontSize',16);
% 
% figure(5)
% plot(t,gtd(1,:),'r-',t,gtd(2,:),'m-',t,gtd(3,:),'b-','linewidth',1);
% hold on
% plot(t,x_esti(1,:),'r--',t,x_esti(2,:),'m--',t,x_esti(3,:),'b--','linewidth',1);
% h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','Location','northwest','NumColumns',3,'FontName','Times New Roman','FontSize',12);
% xlabel('Time','FontName','Times New Roman','FontSize',16);
% ylabel('Position','FontName','Times New Roman','FontSize',16);
% set(h1,'Orientation','horizon','Box','on');
% title('Estimation Result','FontName','Times New Roman','FontSize',16);
% 
% figure(6)
% error = gtd - x_esti;
% error_norm = sqrt(error(1,:).^2 + error(2,:).^2 + error(3,:).^2);
% plot(t,error_norm,'k','linewidth',1);
% xlabel('Time','FontName','Times New Roman','FontSize',16);
% ylabel('Error','FontName','Times New Roman','FontSize',16);
% title('Estimation Error','FontName','Times New Roman','FontSize',16);
% % 
% figure(7)
% % error_x = sqrt(error(1,:).^2 );
% % error_y = sqrt(error(2,:).^2 );
% % error_z = sqrt(error(3,:).^2 );
% plot(t,error(1,:),'r-',t,error(2,:),'m-',t,error(3,:),'b-','linewidth',1);
% xlabel('Time','FontName','Times New Roman','FontSize',16);
% ylabel('Error','FontName','Times New Roman','FontSize',16);
% title('Estimation Error','FontName','Times New Roman','FontSize',16);
% legend('E_{x}','E_{y}','E_{z}','FontName','Times New Roman','FontSize',12);