function val = f_test(x)

A =[...
    2.9770    1.9344   -0.4107
    1.9344    2.8154    0.7750
   -0.4107    0.7750    1.5226];


c = [ ...
    0.1746
   -0.1867
    0.7258];


val = 0.5*x'*A*x + c'*x;









