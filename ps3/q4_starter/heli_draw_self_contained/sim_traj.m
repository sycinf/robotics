function [tsim, xsim] = sim_traj(tout, yout, drawfun, dt)

if(nargin == 3)
    dt = 0.05;
end

minx = min(yout(:,1:3))-.5;
maxx = max(yout(:,1:3))+.5;
axislim = zeros(1,6);
axislim([1 3 5]) = minx;
axislim([2 4 6]) = maxx;

[tsim, xsim] = interp_traj(tout, yout, dt);
for i = 1:length(tsim)
    drawfun(tsim(i), xsim(i,:)', axislim);
end

    function [t_sim, x_sim] = interp_traj(tout, yout, dt)

        t_sim = min(tout):dt:max(tout);
        t_sim = t_sim';
        x_sim = interp1(tout, yout, t_sim);

    end

end