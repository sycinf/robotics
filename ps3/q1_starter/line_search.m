function t = line_search(func,cur_x,cur_dir,alpha,beta,eps)
    t = 1;
    while(func(cur_x+t*cur_dir) > func(cur_x)+alpha*t*gradient_my(func,cur_x,eps).'*cur_dir)
        t=beta*t;
    end
end