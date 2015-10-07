function [x_q_eq_and_ineq_init_min, feasibility] = ...
init_feasible_newton_descent_w_eq_and_ineq...
(fi, x_low, x_high, C, d, x0_feasible, finite_diff_eps, stopping_eps, alpha, beta)
    allFi = fi;
    allFi{length(allFi)+1} = @(x) (x_low - x);
    allFi{length(allFi)+1} = @(x) (x-x_high);
    
    t = 0.5;
    mu = 1.1;
    xinitial = x0_feasible;
    % transform a x+s problem into a y problem
    % we minimize y which is vertcat of x and s
    
    s0_feasible = max(fi{:}(x0_feasible));
    yinitial = vertcat(xinitial,s0_feasible);
    xlen = numel(xinitial);
    coefficientvec = zeros(1, xlen);
    coefficientvec = horzcat(coefficientvec, 1);
    f0=@(y)(coefficientvec*y);
    for consInd = 1:length(fi)
        gcons{consInd}=@(y)(allFi{consInd}(y(1:xlen))-y(xlen+1));
    end
    Cy = horzcat(C,zeros(size(d,1),1))    
    [barrier_y_min, barrier_f_min] = barrier_method...
    (f0,gcons,yinitial,t,mu,Cy,d,finite_diff_eps,stopping_eps,alpha,beta)
    feasibility = (barrier_f_min <= 0);
    x_q_eq_and_ineq_init_min = barrier_y_min(1:xlen);
end