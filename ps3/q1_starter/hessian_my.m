function hes = hessian_my(func,x,eps)
    % how does gradient change when x changes
    centerGrad = gradient_my(func,x,eps);
    gradDim = numel(centerGrad);
    hes = zeros(gradDim,gradDim);
    for cInd = 1: numel(x)
        curXG = x;
        curXL = x;
        curXG(cInd) = curXG(cInd)+eps/2;
        curXL(cInd) = curXL(cInd)-eps/2;
        gradBig = gradient_my(func,curXG,eps);
        gradSmall = gradient_my(func,curXL,eps);
        gradDiff = gradBig-gradSmall;
        hes(:,cInd) = gradDiff/eps;
        
    end
    

end