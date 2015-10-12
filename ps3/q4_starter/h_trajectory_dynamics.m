function [h, jach] = h_trajectory_dynamics(xu_trajectory, cfg)
%function h = h_trajectory_dynamics(xu_trajectory, cfg)
T = cfg.T; nX = cfg.nX; nU = cfg.nU; dt = cfg.dt;

x_traj = reshape(xu_trajectory(1:T*nX),nX, T);
u_traj = reshape(xu_trajectory(T*nX+1:end), nU, T);

f = cfg.f; %dynamics model x_t_plus_1 = f(x_t, u_t, dt);

%% YOUR CODE HERE

%% NOTE: Most equality constraints depend on only a subset of the
%% variables; this makes "blind" numerical computation of the jacobian by
%% the SQP solver very inefficient -- as it'd compute a ton of entries that
%% are known to be zero (but still do all the work for it).  Hence this
%% function also provides the SQP solver with the jacobian of h.  Note that
%% in our experience it is fast enough to compute the non-zero entries in
%% jach through numerical differentation.
h = zeros(nX*T,1);

% here y is vertcat(x,xnext,u)
single_f = @(y)( f(y(1:nX),y(2*nX+1:end),dt)-y(nX+1:2*nX));
jach = zeros(nX*T,T*(nX+nU));

for timeStep = 1:T-1
     curX = x_traj(:,timeStep);
     curU = u_traj(:,timeStep);
     nextX = x_traj(:,timeStep+1);
     
     curY = vertcat(curX,nextX,curU);
     %h(nX*timeStep+1:nX*(timeStep+1)) = single_f(curY);
     h(nX*(timeStep-1)+1:nX*timeStep) = single_f(curY);
     curJac = numerical_jac(single_f,curY);
     %disp('curJac');
     %size(curJac)
%     % now put them inplace
     %jach(nX*timeStep+1:nX*(timeStep+1),(timeStep-1)*nX+1:timeStep*nX) = curJac(:,1:nX);
     %jach(nX*timeStep+1:nX*(timeStep+1),timeStep*nX+1:(timeStep+1)*nX) = curJac(:,nX+1:2*nX);
     %uoffSet = nX*T;
     %jach(nX*timeStep+1:nX*(timeStep+1),uoffSet+(timeStep-1)*nU+1:uoffSet+timeStep*nU)...
     %    = curJac(:,2*nX+1:end);

     
     jach(nX*(timeStep-1)+1:nX*(timeStep),(timeStep-1)*nX+1:timeStep*nX) = curJac(:,1:nX);
     jach(nX*(timeStep-1)+1:nX*(timeStep),timeStep*nX+1:(timeStep+1)*nX) = curJac(:,nX+1:2*nX);
     uoffSet = nX*T;
     jach(nX*(timeStep-1)+1:nX*(timeStep),uoffSet+(timeStep-1)*nU+1:uoffSet+timeStep*nU)...
         = curJac(:,2*nX+1:end);

     
    %f(curX,curU,dt)-nextX;
    %h(nX*timeStep+1:nX*(timeStep+1))= f(x_traj(:,timeStep),u_traj(:,timeStep),dt) - x_traj(:,timeStep+1);
    
end
end