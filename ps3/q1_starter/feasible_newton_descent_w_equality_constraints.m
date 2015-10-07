function [x_min, x_iters, f_iters, stoppingval_iters] = ...
    feasible_newton_descent_w_equality_constraints...
    (f, C, d, x0, finite_diff_eps, stopping_eps, alpha, beta)

% min_x f(x)
% s.t.  Cx = d
%
% newton descent with backtracking line search

% returned values:
% x_min: the minimum found
% x_iters: x values found in each iteration
% f_iters: f values found in each iteration
% stoppingval_iters: delta_x'*H*delta_x from every iteration
    x_iters=[];
    f_iters=[];
    stoppingval_iters=[];
    x_dim = numel(x0);
    cur_x = x0;
    itercount = 0;
    while(true )
        x_iters = vertcat(x_iters,cur_x);
        cur_f = f(cur_x);
        f_iters = vertcat(f_iters,cur_f);
        cur_hes = hessian_my(f,cur_x,finite_diff_eps);
        cur_grad = gradient_my(f,cur_x,finite_diff_eps);
        
        % find the newton step
        LHSCoeMatUpper = horzcat(cur_hes,C');
        vertdimC = size(C,1);
        horzdimC = size(C,2);
        LHSCoeMatLower = horzcat(C,zeros(vertdimC,vertdimC));
        LHSCoeMat = vertcat(LHSCoeMatUpper,LHSCoeMatLower);
        RHS = vertcat(0-cur_grad, zeros(vertdimC ,1));
        solution = inv(LHSCoeMat)*RHS;
        dx_nt = solution(1:x_dim);
        lambda_sq = dx_nt'*inv(cur_hes)*dx_nt;
        stoppingval_iters = vertcat(stoppingval_iters,lambda_sq);
        
        if(lambda_sq/2 <= stopping_eps)
            break;
        end
        t = line_search(f,cur_x,dx_nt,alpha,beta,finite_diff_eps);
        
        cur_x = cur_x+t*dx_nt;
        itercount=itercount+1;
    end
    x_min = cur_x;
    
end
