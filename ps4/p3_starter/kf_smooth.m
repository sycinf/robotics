function [xfilt,varargout] = kf_smooth(y, A, B, C, d, u, Q, R, init_x, init_V)
%
% function [xfilt, xpred, Vfilt, loglik, xsmooth, Vsmooth, Q, R] = 
%           kf_smooth(y, A, B, C, d, u, Q, R, init_x, init_V)
%
%
% Kalman filter
% [xfilt, xpred, Vfilt] = ekf_smooth(y_all, A, B, C, d, Q, R, init_x, init_V);
%
% Kalman filter with Smoother
% [xfilt, xpred, Vfilt, loglik, xsmooth, Vsmooth] = ekf_smooth(y_all, A, B, C, d, Q, R, init_x, init_V);
%
% Kalman filter with Smoother and EM algorithm
% [xfilt, xpred, Vfilt, loglik, xsmooth, Vsmooth, Q, R] = ekf_smooth(y_all, A, B, C, d, Q, R, init_x, init_V);
%
%
% INPUTS:
% y - observations
% A, B, C, d:  x(:,t+1) = A x(:,t) + B u(:,t) + w(:,t) 
%              y(:,t)   = C x(:,t) + d        + v(:,t)
% Q - covariance matrix of system x(t+1)=A*x(t)+w(t) , w(t)~N(0,Q)
% R - covariance matrix of output y(t)=C*x(t)+v(t) , v(t)~N(0,R)
% init_x -
% init_V -
%
%
% OUTPUTS:
% xfilt = E[X_t|t]
% varargout(1) = xpred - the filtered values at time t before measurement
% at time t has been accounted for
% varargout(2) = Vfilt - Cov[X_t|0:t]
% varargout(3) = loglik - loglikelihood
% varargout(4) = xsmooth - E[X_t|0:T]
% varargout(5) = Vsmooth - Cov[X_t|0:T]
% varargout(6) = Q - estimated system covariance according to 1 M step (of EM)
% varargout(7) = R - estimated output covariance according to 1 M step (of EM)


n_var_out = max(nargout,1)-1; % number of variable number of outputs

T = size(y,2);
ss = size(Q,1); % size of state space

%% Forward pass (Filter)
%mu_t_v_0_t = init_x;
%sigma_t_v_0_t = init_V;
%currentu=zeros(size(u,1),1);

% t = 0
loglike = 0;
u0 = zeros(size(u,1),1);
mu_1_v_0_0 = A*init_x+B*u0;
sigma_1_v_0_0 = A*init_V*(A.')+Q;
Vpred(:,:,1) = sigma_1_v_0_0;
K1 = sigma_1_v_0_0*C.'/(C*sigma_1_v_0_0*(C.')+R);
y_0p1 = y(:,1);
mu_1_v_0_1 = mu_1_v_0_0+K1*(y_0p1-(C*mu_1_v_0_0+d));
sigma_1_v_0_1 =  sigma_1_v_0_0- K1*C*sigma_1_v_0_0;

xfilt(:,1)= mu_1_v_0_1;
xpred(:,1)= mu_1_v_0_0;
Vfilt(:,:,1) = sigma_1_v_0_1;

mu_t_v_0_t = mu_1_v_0_1;
sigma_t_v_0_t = sigma_1_v_0_1;
loglike = loglike+curIterLog( sigma_1_v_0_0,...
        mu_1_v_0_0, C, d, R, ss, y_0p1 );

for t = 1:T-1
    % dynamic update
    currentu=u(:,t);
    mu_tp1_v_0_t = A*mu_t_v_0_t + B*currentu;
    sigma_tp1_v_0_t = A*sigma_t_v_0_t*A.'+Q;
    Vpred(:,:,t+1) = sigma_tp1_v_0_t;
    
    % measurement update
    Ktp1 = sigma_tp1_v_0_t*C.'/(C*sigma_tp1_v_0_t*C.'+R);
    y_tp1 = y(:,t+1);
    mu_tp1_v_0_tp1 = mu_tp1_v_0_t+Ktp1*(y_tp1-(C*mu_tp1_v_0_t+d));
    sigma_tp1_v_0_tp1 =  sigma_tp1_v_0_t- Ktp1*C*sigma_tp1_v_0_t;
    
    xfilt(:,t+1)= mu_tp1_v_0_tp1;
    xpred(:,t+1)= mu_tp1_v_0_t;
    Vfilt(:,:,t+1) = sigma_tp1_v_0_tp1;
    
    mu_t_v_0_t= mu_tp1_v_0_tp1;
    sigma_t_v_0_t = sigma_tp1_v_0_tp1;
    loglike = loglike+curIterLog( sigma_tp1_v_0_t,...
        mu_tp1_v_0_t, C, d, R, ss, y_tp1 );
end
disp('xpred');
xpred(:,1)
xpred(:,2)
disp('xfilt');
xfilt(:,1)
xfilt(:,2)
disp('sigma')
Vfilt(:,:,1)
loglike
%YOUR code here

if(n_var_out >= 1), varargout(1) = {xpred}; end
if(n_var_out >= 2), varargout(2) = {Vfilt}; end
if(n_var_out >= 3), varargout(3) = {loglike}; end


%% Backward pass (RTS Smoother and EM algorithm)
xsmooth = zeros(ss,T);
Vsmooth = zeros(ss,ss,T);
if(n_var_out >= 4)
    xsmooth(:,T) = xfilt(:,T);
    Vsmooth(:,:,T) = Vfilt(:,:,T);
    for bt = T-1:-1:1
        Lt = Vfilt(:,:,bt)*A.'*inv(Vpred(:,:,bt+1));
        Vsmooth(:,:,bt) = Vfilt(:,:,bt)+Lt*(Vsmooth(:,:,bt+1)-Vpred(:,:,bt+1))*Lt.';
        xsmooth(:,bt) = xfilt(:,t)+Lt*(xsmooth(:,t+1)-xpred(:,t+1));
    end    
    %YOUR code here
    xsmooth(:,1)
	
	varargout(4) = {xsmooth};
   if(n_var_out >= 5), varargout(4) = {Vsmooth}; end
   if(n_var_out >= 6), varargout(5) = {Q}; end
   if(n_var_out == 7), varargout(6) = {R}; end
end
end


function [logP_y_tp1_v_y_0_t] = curIterLog( sigma_tp1_v_0_t,...
        mu_tp1_v_0_t, C, d, R, dim, y_tp1 )
    y_tp1_v_0_t = C*mu_tp1_v_0_t + d;
    sigma_y_tp1_v_0_t = C*sigma_tp1_v_0_t*C.'+R;
    
    
    lognorm = log(norm(sigma_y_tp1_v_0_t))
    sym =  (y_tp1-y_tp1_v_0_t).'*sigma_y_tp1_v_0_t*(y_tp1-y_tp1_v_0_t)
    
    rtVal = -dim/2*log(2*pi)+1/2*lognorm-1/2*sym;
    
    %p_y_tp1_v_y0t = (2*pi)^(-dim/2) * (norm(sigma_y_tp1_v_0_t))^-0.5 ...
    %    * exp(-0.5*(y_tp1-y_tp1_v_0_t).'*sigma_y_tp1_v_0_t*(y_tp1-y_tp1_v_0_t));
    %p_y_tp1_v_y0t
    
    %-dim/2*log(2*pi)+1/2*log(norm(sigma_y_tp1_v_0_t))
    %logP_y_tp1_v_y_0_t = log(p_y_tp1_v_y0t);
    logP_y_tp1_v_y_0_t = rtVal
end