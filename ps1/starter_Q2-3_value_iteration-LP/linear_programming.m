function [V, pi] = linear_programming(mdp)

% V: optimal value function
% pi: optimal policy
% mdp: description of an MDP
    [numAction] = numel(mdp.T);
    if numAction == 0
        V = zeros(1,1);
        pi = zeros(1,1);
    else
        T = mdp.T;
        R = mdp.R;
        gamma = mdp.gamma;
        array2Dim = T{1};
        [sSize,sPrimeSize] = size(array2Dim);
        mu = ones(sSize,1).';
        cvx_begin
            variable v(sSize)
            minimize (mu*v)
            subject to
            for i=1:sSize 
                for j=1:numAction
                    v(i) >= T{j}(i,:)*(R{j}(i,:)+gamma*v.').'
                end
            end
            V = v;
        cvx_end
        
        % just do one iteration of value iteration
        [lastIterVal, policy] = single_iter_value(mdp,V,numAction,sSize,sPrimeSize);
        pi = policy;
        V = lastIterVal;
        
    end
end