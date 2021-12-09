function [x_esti, x_predict, d] = akf_mhe(gtd, imu, uwb, t)
% This VBAKF-Q use state augmentation method but with a sliding window

%% data preparation
K = length(t);
dt = t(2)-t(1);
x_esti = zeros(8,K);
x_pre = zeros(8,K);
ns = 15; % length of the sliding window
x_esti(1:6,1:ns) = gtd(1:6,1:ns);
x_pre(1:6,1:ns) = gtd(1:6,1:ns);

x_predict = zeros(8,K);
x_predict(1:6,1:ns) = gtd(1:6,1:ns);

y0 = zeros(1,K);
y1 = zeros(1,K);
y0(ns+1:K) = 0.5 * (uwb(1:K-ns).^2);
y1(ns+1:K) = 0.5 * (uwb(ns+1:K).*uwb(ns+1:K));
u = imu;
v = zeros(3,ns+1);
d = zeros(3,K);
y = zeros(1,K);
for i = ns+1:K
    for j = 1:ns+1
    v(1,j) = dt * trapz(u(1,i-ns:i-ns+j-1));
    v(2,j) = dt * trapz(u(2,i-ns:i-ns+j-1));
    v(3,j) = dt * trapz(u(3,i-ns:i-ns+j-1));
    end
    d(1,i) = dt * trapz(v(1,1:ns+1));
    d(2,i) = dt * trapz(v(2,1:ns+1));
    d(3,i) = dt * trapz(v(3,1:ns+1));
    y(i) = y1(i) - y0(i) + 0.5 * (d(1,i)^2 + d(2,i)^2 + d(3,i)^2);
end

%% x_k = A x_k-1 + B u_k
% Should A and B change when using MHE?
A = [1,0,0,dt,0,0,0,0;
     0,1,0,0,dt,0,0,0;
     0,0,1,0,0,dt,0,0;
     0,0,0, 1,0,0,0,0;
     0,0,0, 0,1,0,0,0;
     0,0,0, 0,0,1,0,0;
     0,0,0, 0,0,0,1,0;
     0,0,0, 0,0,0,0,1];

% p0,vo is changing
A1 = [1,0,0,dt,0,0,0,0;
     0,1,0,0,dt,0,0,0;
     0,0,1,0,0,dt,0,0;
     0,0,0, 1,0,0,0,0;
     0,0,0, 0,1,0,0,0;
     0,0,0, 0,0,1,0,0;
     0,0,0, 0,0,0,0,0;
     0,0,0, 0,0,0,0,0];

B = [0.5*dt^2,0,0;
     0,0.5*dt^2,0;
     0,0,0.5*dt^2;
     dt,   0,   0;
     0,   dt,   0;
     0,    0,  dt;
     0,    0,   0;
     0,    0,   0];
 
%% Parameter
rho = 1 - 1e-4;
nx = 8;
tau  = 2;
R = 1e-2;
P = diag([1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-2,1e-1]);
Q0 = diag([1e-1,1e-1,1e-1,1e-7,1e-7,1e-7,1e-2,1e-4]);
mu = tau + nx +1;
U = tau * Q0;
delta = 1e-6;
N = 10;

%% VBAKF-Q + MHE + State Augmentation
for i = ns+1:K
%     disp(['Step: ',int2str(i)]);
    x_predict(:,i) = A1 * x_predict(:,i-1) + B * u(:,i);
    x_predict(7,i) = x_predict(1:3,i-ns)'*x_predict(4:6,i-ns);
    x_predict(8,i) = x_predict(4:6,i-ns)'*x_predict(4:6,i-ns);
    % Prediction
    x_pre(:,i) = A1 * x_esti(:,i-1) + B * u(:,i);
    x_pre(7,i) = x_esti(1:3,i-ns)'*x_esti(4:6,i-ns);
    x_pre(8,i) = x_esti(4:6,i-ns)'*x_esti(4:6,i-ns);
    sigma = A * P * A';
    mu = rho * (mu-nx-1) + nx + 1;
    U = rho * U;
    % Update
    H = [d(1,i),d(2,i),d(3,i),0,0,0,(ns * dt),0.5*((ns * dt)^2)];
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
    % if NaN in x_esti?
    TF = isnan(x_esti(:,i));
    if ismember(1,TF)
        disp(['Divergence! The step is: ',int2str(i)]);
        break;
    end
    
end

