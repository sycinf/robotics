function [x_q_g_min, x_q_g_iters, f_q_g_iters, stoppingvals_q_g_iters] = ...
    gradient_descent(func, x, finite_diff_eps, stopping_eps, alpha, beta)
    cur_grad = gradient_my(func,x,finite_diff_eps);
    cur_x = x;
    cur_f=func(cur_x);
    x_q_g_iters=cur_x;
    f_q_g_iters=cur_f;
    stoppingvals_q_g_iters=[];
    while(norm(cur_grad,2) >= stopping_eps)
        cur_dir = 0-cur_grad;
        t = line_search(func,cur_x,cur_dir,alpha,beta,finite_diff_eps);
        cur_x = cur_x + cur_dir*t;
        cur_grad = gradient_my(func,cur_x,finite_diff_eps);
        cur_f=func(cur_x);
        x_q_g_iters = vertcat(x_q_g_iters,cur_x);
        f_q_g_iters = vertcat(f_q_g_iters,cur_f);
    end
    x_q_g_min = cur_x;
    
end