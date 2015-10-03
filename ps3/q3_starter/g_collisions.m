function [val,jac] = g_collisions(x, dsafe, traj_shape, make_robot_poly, obstacles)
assert(size(x,2)==1);
traj = reshape(x, traj_shape);
[K,T] = size(traj);

val = zeros(length(obstacles)*size(traj,2),1);
jac = zeros(size(val,1), size(x,1));

icontact = 1;


for t=1:T
    xt = traj(:,t);
    for iobs=1:length(obstacles)
        [d,pts] = signedDistancePolygons(...
                make_robot_poly(xt), ...
                obstacles{iobs});
        ptOnRobot = pts(1,:);
        ptOnObs = pts(2,:);
        normalObsToRobot = -sign(d)*normr(ptOnRobot - ptOnObs);

        gradd = normalObsToRobot * calcJacobian(ptOnRobot, xt);
        
        val(icontact) = dsafe - d;
        jac(icontact,K*(t-1)+1:K*t) = gradd;
        icontact = icontact+1; 
   end 
end


end

function jac = calcJacobian(pt, x0)
    jac = zeros(2,3);
    r = pt(:) - x0(1:2);
    jac(1,1) = 1;
    jac(2,2) = 1;
    jac(1,3) = -r(2);
    jac(2,3) = r(1);
end