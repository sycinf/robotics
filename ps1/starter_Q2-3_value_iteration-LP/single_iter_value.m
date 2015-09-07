function [perIterV,actionInd] = single_iter_value(mdp,prevV,numAction,sSize,sPrimeSize)

    curIterV = zeros(sSize,numAction);
    
    for aInd = 1:numAction
        for sInd = 1:sSize
            for sPInd = 1:sPrimeSize
                if  mdp.T{aInd}(sInd,sPInd) > 0
                    curIterV(sInd,aInd) = curIterV(sInd,aInd) + mdp.T{aInd}(sInd,sPInd)*...
                    (mdp.R{aInd}(sInd,sPInd)+mdp.gamma*prevV(sPInd));
                end
            end
        end
    end
    % traverse the curV array to get the max
    [perIterV,actionInd] = max(curIterV.');
    perIterV = perIterV.';
    actionInd = actionInd.';
end
