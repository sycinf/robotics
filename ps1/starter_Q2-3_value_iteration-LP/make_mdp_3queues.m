function [ mdp ] = make_mdp_3queues(addReqProb,reward, qSize, gamma)
    mdp.gamma = gamma;
    totalNumStates = (qSize+1)^3*3
    T{1} = zeros(totalNumStates,totalNumStates);
    T{2} = zeros(totalNumStates,totalNumStates);
    T{3} = zeros(totalNumStates,totalNumStates);
    
    R{1} = zeros(totalNumStates,totalNumStates);
    R{2} = zeros(totalNumStates,totalNumStates);
    R{3} = zeros(totalNumStates,totalNumStates);
    % a lot of state to generate
    maxStateNumber = 0;
    minStateNumber = 500000;
    for actionInd = 1:3
        for aQueueSize = 0:qSize
            for bQueueSize = 0:qSize
                for cQueueSize = 0:qSize
                    for queuePosition = 0:2
                        curState = queueState;
                        curState.queueSize = zeros(3,1);
                        curState.queueSize(1) = aQueueSize;
                        curState.queueSize(2) = bQueueSize;
                        curState.queueSize(3) = cQueueSize;
                        curState.queuePos = queuePosition;
                        curStateNumber = curState.getStateNumeric(qSize);
                        [obtainedReward,allNextState, allNextStateProb] = curState.generateNextStateT(addReqProb,qSize,actionInd, reward);
                        numNextState = numel(allNextState);
                        for nextStateInd = 1:numNextState
                            curNextState = allNextState(nextStateInd);
                            
                            curNextStateProb = allNextStateProb(nextStateInd);
                            if curNextStateProb ~= 0
                                curNextStateNumber = curNextState.getStateNumeric(qSize);
                                T{actionInd}(curStateNumber,curNextStateNumber) = curNextStateProb;
                               
                                if obtainedReward
                                    R{actionInd}(curStateNumber,curNextStateNumber) = obtainedReward;
                                end
                            end
                        end  
                        if curStateNumber >maxStateNumber
                            maxStateNumber = curStateNumber;
                        end
                        if curStateNumber < minStateNumber
                            minStateNumber = curStateNumber;
                        end
                    end
                end
            end
        end
    end
    
    mdp.T = T;
    mdp.R = R;
   
        
        
        
end

