classdef (Abstract) STDPModel_delays < STDPModel
  %SynapseModel_g_exp Conductance-based single exponential synapses
  %   Parameters to set in ConnectionParams:
  %   - E_reversal, the reversal potential (in mV)
  %   - tau, the synaptic decay time constant (in ms)


    
  properties (SetAccess=private)
    STDPbufferCount
    STDPbufferMax
    preBoundaryArrPerN
    preABSID
  end
  
  properties (SetAccess=public)
    preGroups
    spikeArrLogicInd
  end
  
  methods
      % Synapse model will be a cell array containing a struct for each
      % connection type indexed by the post synaptic group then the
      % presynaptic synapse group(may contain multiple neuron groups). The struct contains variables for all
      % neurons in the presynatpic synapse group and all neurons in the post
      % synaptic group. Or all post synaptic neurons on the lab and all
      % presynaptic neurons regardless of lab, if using parallel mode. 
    function SM = STDPModel_delays(CP, SimulationSettings, ...
                                     postID, number_in_post,number_in_pre,pre_group_ids, group_boundary)
    
      SM = SM@STDPModel(CP, SimulationSettings, ...
                                     postID, number_in_post,number_in_pre,pre_group_ids);


      SM.STDPbufferCount = 1;
      
      SM.preBoundaryArrPerN = zeros(sum(number_in_pre),2,'uint16');
      SM.preABSID = zeros(sum(number_in_pre),1, 'uint16');

      for iB = 2:length(SM.preBoundaryArr)
          
          SM.preBoundaryArrPerN(SM.preBoundaryArr(iB-1)+1:SM.preBoundaryArr(iB),1 ) = SM.preBoundaryArr(iB-1);
          SM.preBoundaryArrPerN(SM.preBoundaryArr(iB-1)+1:SM.preBoundaryArr(iB),2 ) = SM.preGroupIDs(iB-1);
          SM.preABSID(SM.preBoundaryArr(iB-1)+1:SM.preBoundaryArr(iB),1 ) = ...
              (SM.preBoundaryArr(iB-1)+1:SM.preBoundaryArr(iB))-SM.preBoundaryArr(iB-1) ...
              +group_boundary(SM.preGroupIDs(iB-1));
      end
      
      
      maxDelaySteps = SimulationSettings.maxDelaySteps;
      SM.spikeArrLogicInd = cell(maxDelaySteps,1);
      
      for i = 1:length(SM.spikeArrLogicInd)
          SM.spikeArrLogicInd{i} = speye(number_in_post,sum(number_in_pre));
          SM.spikeArrLogicInd{i}(:) = false;
      end
      %trace variable for presynaptic neurons, contains an entry for each
      %neuron in the presynaptic group of the connection.
      SM.Apre = zeros(maxDelaySteps,sum(number_in_pre));

      
      SM.STDPbufferMax = maxDelaySteps;

    end

    
    function SM = updateSTDPBuffer(SM)
      
      %SM.Apre(SM.bufferCount,:) = 0;
      SM.Apre(SM.STDPbufferCount+1,:) = SM.Apre(SM.STDPbufferCount,:);
      SM.spikeArrLogicInd{SM.STDPbufferCount}(:) = false;
      
      SM.STDPbufferCount = SM.STDPbufferCount + 1;
      
      if SM.STDPbufferCount > SM.STDPbufferMax
        SM.STDPbufferCount = 1;
        SM.Apre(SM.STDPbufferCount,:) = SM.Apre(SM.STDPbufferMax,:);
      end
    end
    
    function SM = updateSynapses(SM, dt)
        

      SM.Apre(SM.STDPbufferCount,:) = SM.Apre(SM.STDPbufferCount,:) + (( - SM.Apre(SM.STDPbufferCount,:))./SM.tPre).*dt;
      SM.Apost = SM.Apost + ((-  SM.Apost)./SM.tPost).*dt;
      
    end
    
    function weightsArr = updateWeights(SM,weightsArr,IDMap,postGroup, synArr)
        %Retreive relative pre and post synaptic IDs of arriving spikes
        [post, pre] = find(SM.spikeArrLogicInd{SM.STDPbufferCount});
        %Convert to absolute IDs for accessing weights array.
        postAbs = IDMap.cellIDToModelIDMap{postGroup}(post);
        pre = SM.preABSID(pre);
        for i = find(abs(SM.Apost(post))>0)
            postLoc = (synArr{pre(i), 1}== postAbs(i));
            weightsArr{pre(i)}(postLoc) = weightsArr{pre(i)}(postLoc) + SM.Apost(post(i));
        end
    end

    

    
    % Called when spike has been generated by a post synaptic neuron in this group 
    function weightsArr = updateweightsaspostsynspike(SM, weightsArr, preInd,groupID, delay)
        %update the weight array for all connections from spiking neuron to
        %all pre synaptic neurons.

        preInd = preInd + SM.preBoundaryArr(ismember(SM.preGroupIDs,groupID));

        ind = sub2ind(size(SM.Apre), delay, preInd);


        weightsArr = weightsArr + SM.Apre(ind)';

        weightsArr(weightsArr<SM.wmin) = SM.wmin;
        weightsArr(weightsArr>SM.wmax) = SM.wmax;
    end
    
    %spikeInd are the indices of presynaptic neurons fromt this group
    %that have spiked during this cycle
    function SM = processAsPreSynSpike(SM, preSpikeInd,groupID, postInd,bufferLoc)
        %update presynaptic trace variable Apre, should be same for all
        %presynaptic connections because the synapse parameters should be
        %the same for all group to group defined connections. 
        
        preSpikeInd = preSpikeInd + SM.preBoundaryArr(SM.preGroupIDs==groupID); % convert from neuron group relative ID to synapse group relative ID

        %for each post synaptic neuron
        for i = 1:length(bufferLoc)
            SM.spikeArrLogicInd{bufferLoc(i)}(postInd(i),preSpikeInd) = true;
        end% place into the spike buffer at the buffer location (the current time + delay)
        
        SM.Apre(SM.STDPbufferCount,preSpikeInd) = SM.Apre(SM.STDPbufferCount,preSpikeInd) + SM.preRate; % Increment Apre at the current position. 
    end
    


    

    
    function check = hasArrivingSpikes(SM)
        check = nnz(SM.spikeArrLogicInd{SM.bufferCount})>1;
    end
    

    

    function [Apre] = getApre(SM, preInd, groupID)

        preInd = preInd + SM.preBoundaryArr(SM.preGroupIDs==groupID);
        Apre = SM.Apre(SM.STDPbufferCount,preInd);
    end

  end % methods
  

end % classdef
