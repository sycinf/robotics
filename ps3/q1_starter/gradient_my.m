function grad = gradient_my(func,x, eps)
    
    sizeGradC = numel(x);
    sizeGradR = numel(func(x));
    grad = zeros(sizeGradR, sizeGradC);
    curXG = x;
    curXL = x;
    % generating x one by one
    for cInd = 1:sizeGradC
        curXG(cInd)=curXG(cInd)+eps/2;
        curXL(cInd)=curXL(cInd)-eps/2;
        deltaF = func(curXG) - func(curXL);
        grad(:,cInd) = deltaF/eps;
        curXG = x;
        curXL = x;
    end        
    grad=grad';
end