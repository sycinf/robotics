function [barrier_x_min, barrier_f_min] = barrier_method...
    (f0,fcons,x0,t,mu,C,d,finite_diff_eps,stopping_eps,alpha,beta)
    % x0 is strictly feasible, t > 0, mu >1, stopping_eps > 0 
    cur_t = t;
    cur_x = x0;
    while(true)
        centeringFunc = @(x)(newFunc(f0,cur_t,fcons,x));
        [x_min, x_iters, f_iters, stoppingval_iters] = ...
            feasible_newton_descent_w_equality_constraints...
            (centeringFunc, C, d, cur_x, finite_diff_eps, stopping_eps, alpha, beta);
        cur_x=x_min;
        if length(fcons)/cur_t < stopping_eps
            break;
        end
        cur_t = mu*t;
    end
    barrier_x_min = cur_x;
    barrier_f_min = f0(cur_x);
end

function fresult = newFunc(f0,cur_t,fcons,x)
    %numCons = length(fcons);
    %sumLogCon = 0;
    %for consInd = 1: numCons
    %    curFcon = fcons{consInd};
    %    sumLogCon = sumLogCon+ log(-curFcon(x));
    %end
    %fresult = f0(x)-(1/t)*sumLogCon;
    fresult=f0(x)-(1/cur_t)* log((-fcons{:}(x))')*ones(length(fcons),1);
end
