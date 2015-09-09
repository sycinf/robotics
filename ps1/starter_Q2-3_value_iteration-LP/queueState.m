classdef queueState
    % state space a0...a5 x b0...b5 x c0....c5
    % 3*216 state
    % a0...a5 would be 0....5 
    % b0...b5 would be 0....5 * 6
    % c0...c5 would be 0....5 * 36
    % cur position in a -- prev + 0*216
    % cur position in b -- prev + 1*216
    % cur position in c -- prev + 2*216
    % add a offset one
    % queueSize is a vector of 3
    properties
        queueSize
        queuePos        
    end
    
    methods
        function rtNum = getStateNumeric(obj,maxQueueSize)
            rtNum = (obj.queueSize(1)+...
                obj.queueSize(2)*(maxQueueSize+1)+...
                obj.queueSize(3)*(maxQueueSize+1)^2)+...
                obj.queuePos*(maxQueueSize+1)^3 + 1;
        end
        % action = 1 --> stay, action = 2 --> move left, action = 3 -->
        % move right
        function [obtainedReward,nextState, nextStateProb] = generateNextStateT(obj,addReqProb,maxQueueSize,action, reward)
            queueSizeB4Add = obj.queueSize;
            if action == 1 && queueSizeB4Add(obj.queuePos+1) > 0
                queueSizeB4Add(obj.queuePos+1) = queueSizeB4Add(obj.queuePos+1) - 1;
                obtainedReward = reward(obj.queuePos+1);
            else
                obtainedReward = 0;
            end

                
            % from here, we can generate 8 different state    
            % some of which are non-valid, so we will traverse them
            % after they are created, and then for those impossible ones
            % we take their probability and give it to the ones which will
            % be possible -- any state with queue size goes up 6, maps to
            % the ones with all the queues' size = 5
            nextStateProb = zeros(8,1);
            for aAdd = 0:1
                for bAdd = 0:1
                    for cAdd = 0:1
                        stateIndex = aAdd+bAdd*2+cAdd*4+1;
                        curProbability = gen_probability(addReqProb,[aAdd,bAdd,cAdd]);
                        curNextState = queueState;
                        curNextState.queueSize = queueSizeB4Add;
                        
                        curNextState.queueSize(1) = curNextState.queueSize(1)+aAdd;
                        curNextState.queueSize(2) = curNextState.queueSize(2)+bAdd;
                        curNextState.queueSize(3) = curNextState.queueSize(3)+cAdd;
                        curNextState.queuePos = obj.queuePos;
                        if action == 2
                            curNextState.queuePos = curNextState.queuePos - 1;
                            if curNextState.queuePos == -1
                                curNextState.queuePos = 2;
                            end
                        end
                        if action == 3
                            curNextState.queuePos = curNextState.queuePos + 1;
                            if curNextState.queuePos == 3
                                curNextState.queuePos = 0;
                            end
                        end
                        
                        nextState(stateIndex) = curNextState;
                        nextStateProb(stateIndex) = curProbability;
                    end
                end
            end
            % scan all the next state to see if they are legit
            for stateInd = 1:8
                curNextState = nextState(stateInd);
                curNextStateProb = nextStateProb(stateInd);
                complementaryStateInd = stateInd;
                for queueInd = 1:3
                    if curNextState.queueSize(queueInd) > maxQueueSize
                        complementaryStateInd = complementaryStateInd - 2^(queueInd-1);
                    end
                end
                % not a legit state, give all their probability to the
                %
                if complementaryStateInd ~= stateInd
                    nextStateProb(stateInd) = 0;
                    nextStateProb(complementaryStateInd) = ...
                        nextStateProb(complementaryStateInd)+curNextStateProb;
                end                
            end 
        end
        
        
    end
end

function probability = gen_probability(addReqProb,queueIncremented)
    probability = 1;
    for queueInd = 1:3
        curAddReqProb = addReqProb(queueInd);
        if queueIncremented(queueInd) == 0
            probability = probability*(1-curAddReqProb);
        else
            probability = probability*curAddReqProb;
        end           
    end
end
