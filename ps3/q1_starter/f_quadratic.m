function val = f_quadratic(x)

c = [0.015273927029036; 0.746785676564429; 0.445096432287947];

A = [     2.976991824249534   1.934366257364652  -0.410688884621065
   1.934366257364652   2.815433927071369   0.775017509184702
  -0.410688884621065   0.775017509184702   1.522648051834486];

%grad:A*x+c
%hessian: A
val = 0.5*x'*A*x + c'*x;