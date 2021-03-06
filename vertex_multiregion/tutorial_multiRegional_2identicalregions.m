
% set up the two identical layers by calling setup_multilayer and cloning
% the parameters for the second region.

 setup_multiregion_withinboundconnection;

RecordingSettings.saveDir = '~/VERTEX_results_2regionconnected/stim1000R1toR2';
RecordingSettings.weights_arr = 1:5000:100000;
SimulationSettings.simulationTime = 500;
SimulationSettings.timeStep = 0.03125;
SimulationSettings.parallelSim = true;
% 
[TissueParams.StimulationField,model] = invitroSliceStim('tutorial2_3.stl',100);
TissueParams.StimulationOn = 300;
TissueParams.StimulationOff = 350;%SimulationSettings.simulationTime;

% optional - step current stimulation to neurons to see spread of activity
% through region to region connections
% 
%  NeuronParams(1).Input(2).inputType = 'i_step';
%  NeuronParams(1).Input(2).timeOn = 50;
%  NeuronParams(1).Input(2).timeOff = 100;
%  NeuronParams(1).Input(2).amplitude = 1000; 

[params, connections, electrodes] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams, ...
              RecordingSettings, SimulationSettings);

 % clone the slice to create an identical second region.\
 % (if wanting two differing regions you will need to call a second version
 % of setup_multilayer with mofified parameters and then initialise this
 % new network with another call of initNetwork. As commented out below:
 
 clear TissueParams NeuronParams ConnectionParams
 
 
 setup_multiregion_withinboundconnection;

 [params2, connections2, electrodes2] = ...
  initNetwork(TissueParams, NeuronParams, ConnectionParams, ...
              RecordingSettings, SimulationSettings);
 

 

%%  
 % defining the between region connectivity here. If you have connection (fibre) lengths
 % in mm then a matrix of fibre lengths can also be passed here, and delays will
 % be calculated automatically, assuming a transmission speed of 120mm/ms
 
 regionConnect.map = [0,5000;0,0];
 % for example [1,1;0,1] there are two regions and there is only an
 % external connection from region 1 to region 2, it is not returned, and
 % while they do connect to themselves internally for the sake of incoming external
 % connections the diagonals are set to 0.
 
 % Identify the neuron types (e.g. NP(1) in this instance are the
 % excitatory PY cells) which will export signals. Use [] if not exporting.
 regionConnect.exportingNeuronPops{1} = 1; 
 regionConnect.exportingNeuronPops{2} = 1;
 % identify which neuron pops are designated as dummy neurons to just
 % recieve external signals. (May need to change this if having a different
 % implementation of the dummy neurons.) Use [] if no dummy neurons are
 % present.
 regionConnect.dummyNeuronPops{1} = [];
 regionConnect.dummyNeuronPops{2} = 3;

 
 %stack the parameters for params, connections and electrodes into cell
 %arrrays which can then be fed into runSimulationMultiregional
 paramStacked = {params, params2};
 connectionStacked = {connections,connections2};
 electrodeStacked = {electrodes,electrodes2};
 
runSimulationMultiregional(paramStacked,connectionStacked,electrodeStacked,regionConnect);

% need to use a Multiregions variant of loadResults to load the results for
% every region in one structure. 
Results = loadResultsMultiregions(RecordingSettings.saveDir);
% 
 
[TissueParams.StimulationField,model] = invitroSliceStim('tutorial2_3.stl',100);
figure; 
pdeplot3D(model,'ColorMapData', TissueParams.StimulationField.NodalSolution, 'FaceAlpha', 0.2);
hold on
plotSomaPositions(Results(1).params.TissueParams)
hold on
scatter3(RecordingSettings.meaXpositions(:),RecordingSettings.meaYpositions(:),RecordingSettings.meaZpositions(:),'*')

 
plotSpikeRaster(Results(1))
title('Region 1')
plotSpikeRaster(Results(2))
title('Region 2')
figure
plot(mean(Results(1).LFP,2))
hold on
plot(mean(Results(2).LFP,2))

