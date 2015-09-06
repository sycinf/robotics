% value iteration:
function [V, pi] = value_iteration(mdp, precision)

%IN: mdp, precision
%OUT: V, pi

% Recall: to obtain an estimate of the value function within accuracy of
% "precision" it suffices that one of the following conditions is met:
%   (i)  max(abs(V_new-V)) <= precision / (2*gamma/(1-gamma))
%   (ii) gamma^i * Rmax / (1-gamma) <= precision  -- with i the value
%   iteration count, and Rmax = max_{s,a,s'} | R(s,a,s') |
    [numAction] = numel(mdp.T);
    if numAction == 0
            V = zeros(1,1);
            pi = zeros(1,1);
    else
            array2Dim = (mdp.T{1});
            [sSize,sPrimeSize] = size(array2Dim);
            prevV = zeros(sSize,1);
            curV = single_iter_value(mdp,prevV,numAction,sSize,sPrimeSize);
            gamma = mdp.gamma;
            
            Rmax = max(max(max(cat(3, mdp.R{:}))))
            i = 0;
            while (max(abs(curV-prevV))>(precision / (2*gamma/(1-gamma))))...
                    && (gamma^i * Rmax/(1-gamma) > precision )
                prevV = curV;
                [curV, pi] = single_iter_value(mdp,prevV,numAction,sSize,sPrimeSize);
            end
            V = curV
            
    end
end

function [perIterV,actionInd] = single_iter_value(mdp,prevV,numAction,sSize,sPrimeSize)

    curIterV = zeros(sSize,numAction);
    disp(curIterV(1))
    for aInd = 1:numAction
        for sInd = 1:sSize
            for sPInd = 1:sPrimeSize
                if  mdp.T{aInd}(sInd,sPInd) > 0
                    curIterV(sInd,aInd) = curIterV(sInd,aInd) + mdp.T{aInd}(sInd,sPInd)*...
                    (mdp.R{aInd}(sInd,sPInd)+mdp.gamma*prevV(sPInd))
                end
            end
        end
    end
    % traverse the curV array to get the max
    [perIterV,actionInd] = max(curIterV.');
    perIterV = perIterV.';
    actionInd = actionInd.';
end


