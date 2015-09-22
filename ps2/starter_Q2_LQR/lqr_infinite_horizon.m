function [K_inf, P_inf] = lqr_infinite_horizon(A, B, Q, R)
    xDim = size(A,1)
    P0 = zeros(xDim,xDim)
    Pprev = P0;
    
    K1 = -inv(R+B.'*Pprev*B)*B.'*Pprev*A;
    Kprev = K1;
    Kcur = K1;
    
    Pcur = Q+Kcur.'*R*Kcur+(A+B*Kcur).'*Pprev*(A+B*Kcur); %P1
    while(1)
        Kprev = Kcur;
        Pprev = Pcur;
        Kcur =  -inv(R+B.'*Pprev*B)*B.'*Pprev*A;
        Pcur = Q+Kcur.'*R*Kcur+(A+B*Kcur).'*Pprev*(A+B*Kcur);             
        if (norm(Kcur - Kprev, 2) <= 1e-4)
            break;
        end
    end
    
    K_inf = Kcur
    P_inf = Pcur
end
