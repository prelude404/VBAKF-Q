function [x_esti,x_predict] = vbakf_q(gtd, imu, uwb, t)
% vbakf_q: The variational Bayesian adaptive Kalman filter with unknown Q.
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
x_predict = zeros(8,K);
x_esti(:,1) = [x0(1:3,1); x0(4:6,1); x0(1:3,1)'*x0(4:6,1); x0(4:6,1)'*x0(4:6,1)];
x_predict(:,1) = x_esti(:,1);
x_pre = zeros(8,K);
x_pre(:,1) = x_esti(:,1);

%% Parameter
rho = 1 - 1e-4;
nx = 8;
tau = 2;
% P = diag([1e-2*[1,1,1],1e-3*[1,1,1],1e-2,1e-1]);
% Q0 = diag([1e-1*[1,1,1],1e-7*[1,1,1],1e-2,1e-4]);
R = 1e-10;
% Q0 = 1e-5 * [0.0001*eye(3),0.0032*eye(3),zeros(3,2);
%             0.0032*eye(3),0.1600*eye(3),zeros(3,2);
%             zeros(2,8)];

P = diag([1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-2,1e-1]);
Q0 = diag([1e-1,1e-1,1e-1,1e-7,1e-7,1e-7,1e-2,1e-4]);


% Q0 = diag([1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-6].^2);
% R = 100;
% P = diag([10*ones(1,3) 0.01*ones(1,3) 1 1])*0.001;

%% Parameter for data bag in real flight
% rho = 1 - 1e-4;
% nx = 8;
% tau = 2;
% P = diag([1e1*[1,1,1],1e-3*[1,1,1],1e-2,1e-5]);
% Q0 = diag([1e1*[1,1,1],1e-7*[1,1,1],1e-4,1e-6])*0.1;
% R = 1e2;

mu = tau + nx +1;
U = tau * Q0;
delta = 1e-6;
N = 10;

%% VBAKF-Q
for i = 2:K
%     disp(['Step: ',int2str(i)]);
    % Prediction
    x_predict(:,i) = A * x_predict(:,i-1) + B * u(:,i);
    x_pre(:,i) = A * x_esti(:,i-1) + B * u(:,i);
    sigma = A * P * A';
    mu = rho * (mu-nx-1) + nx + 1;
    U = rho * U;
    % Update
    H = [d(1,i),d(2,i),d(3,i),0,0,0,((i-1)*dt),0.5*(((i-1)*dt)^2)];
    theta = x_pre(:,i);
    x_former = theta;
    for j = 1:N
        Aq = U./mu;
        % update x
        K = Aq * H' * ((H*Aq*H'+R)^(-1));
        x = theta + K * (y(:,i) - H * theta);
        P = Aq - K * H * Aq;
        % update theta
        K1 = sigma * ((sigma + Aq)^(-1));
        theta = x_pre(:,i) + K1 * (x - x_pre(:,i));
        P1 = sigma - K1 * sigma;
        % update Q
        mu = mu + 1;
        U = U + (x-theta)*(x-theta)' + P + P1;
        if norm(x - x_former)/norm(x) < delta
            break;
        end
        x_former = x;
    end
    x_esti(:,i) = x;
    
%     % if NaN in x_esti?
%     TF = isnan(x_esti(:,i));
%     if ismember(1,TF)
%         disp(['Divergence! The step is: ',int2str(i)]);
%         break;
%     end
    
end

