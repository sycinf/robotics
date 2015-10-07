function [x_q_newton_min, x_q_newton_iters, f_q_newton_iters, stoppingvals_q_newton_iters]...
    = newton_descent( func, x, finite_diff_eps, stopping_eps, alpha, beta)
    x_q_newton_iters=[];
    f_q_newton_iters=[];
    stoppingvals_q_newton_iters=[];
    cur_x = x;
    while(true)
        x_q_newton_iters = vertcat(x_q_newton_iters,cur_x);
        cur_f = func(cur_x);
        f_q_newton_iters = vertcat(f_q_newton_iters,cur_f);
        cur_hes = hessian_my(func,cur_x,finite_diff_eps);
        cur_grad = gradient_my(func,cur_x,finite_diff_eps);
        dx_nt = 0-inv(cur_hes)*cur_grad;
        %dx_nt= 0-cur_hes\cur_grad;
        lambda_sq = cur_grad'*inv(cur_hes)*cur_grad;
        %lambda_sq = cur_grad'/cur_hes*cur_grad;
        
        if(lambda_sq/2 <= stopping_eps)
            break;
        end
        t = line_search(func,cur_x,dx_nt,alpha,beta,finite_diff_eps);
        
        cur_x = cur_x+t*dx_nt;
    end
    x_q_newton_min = cur_x;
    
end