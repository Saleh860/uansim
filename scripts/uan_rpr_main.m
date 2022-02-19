function [results,t]=uan_rpr_main(basedir)
	if nargin<1
		basedir='..\';
	end

	dispTopology=false;
	debugSimulation=0;

	javapath=[basedir,'bin'];
    dpath=javaclasspath();
    if sum(strcmpi(dpath,javapath))>0
        javarmpath(javapath);
    end
    javaaddpath(javapath);
    
	octavepath=[basedir,'scripts\octave'];
    addpath(octavepath);

	% load("pos_data.mat", "pos_data");
	% pos_data=pos_data(1,:);
	% pos_data=pos_data(:,1:10);

	load([basedir, 'data\good_pos_data.mat'], "good_pos_data");
	pos_data=good_pos_data(:,1:100);

	%   network_density number_of_nodes min_number_of_neighbors
	experiments = [
		4	28	3
		5	47	4
	...    6	57	5
	...    7	66	6
		8	75	7
	...    9	85	8
		15	140	9
	];

	%Max # of attackers
	A=1;

	%Result columns
	z=zeros(size(pos_data,1)*(1+A),1);
	D=z; p=z; attacks=z;nodes=z; src_deg=z;genrtd=z;sent=z;recved=z;forwd=z;pktcost=z;
	e2edly=z; maxpktcost=z; fintim=z;

	r=0;
	% for each experiment
	for ii = 1: size(pos_data,1) 
		%for each number of attackers
		for a  = 0:A  
			r=r+1;
			D(r)=experiments(mod(ii-1, length(experiments))+1,1);
			attacks(r)=a;
			
			%number of topologies
			n=size(pos_data,2);
			z=zeros(n,1);
			s_nodes=z;
			s_src_deg=z;
			s_genrtd=z;
			s_sent=z; 
			s_recved=z;
			s_forwd=z;
			s_pktcost=z;
			s_e2edly=z;
			s_maxpktcost=z;
			s_fintim=z;
			
			for k= 1:n
				sim=createJavaObject('uansim.RPRSimulator');

				%%%%%%%%%% PARAMETERS %%%%%%%%%%
				% Topology parameters
				sim.topologyParameters.topologySeed = 28629;
				sim.topologyParameters.networkDensity = experiments(mod(ii-1, size(experiments,1))+1,1);
				sim.topologyParameters.deploymentSideLength=500;
				sim.topologyParameters.deploymentDepth=250;
				sim.topologyParameters.sourcesCount = 1;
				sim.topologyParameters.attackersCount = a;
				nodeCount = experiments(mod(ii-1, length(experiments))+1,2);

				% Traffic parameters
				sim.trafficParameters.messageCount=100;
				sim.trafficParameters.generationBitRate=1;
				sim.trafficParameters.messageBitLength=50*8;
				sim.trafficParameters.channelBitRate=100000/8;

				% Routing parameters
				sim.routingParameters.K = .0;
				sim.routingParameters.forwardingThreshold = 0;
				sim.routingParameters.maximumForwardingDelay=1;
				sim.routingParameters.decisionSeed = 64641;
				sim.routingParameters.discoveryRounds = 5;
				sim.routingParameters.maximumTTL = 20;

				%Simulation parameters
				sim.debug = 0;
				dispTopology = false;
				%%%%%%%%%% Create nodes %%%%%%%%%%%
				%sim.topology.prepare(); 
				createTopology(sim, pos_data{ii,k});
				if dispTopology
					disp('sim.topology.sinkPos=new double[]')
					t=sim.topology.sinkPos;
					disp(['{',mat2str(t(1)),',',mat2str(t(2)),',',mat2str(t(3)),'};']);
					disp('sim.topology.nodePositions=new double[][]{')
					t=sim.topology.nodePositions;
					for i=1:size(t,1)
						disp(['{',mat2str(t(i,1)),',',mat2str(t(i,2)),',',mat2str(t(i,3)),'},']);
					end
					disp('};')
					disp('sim.topology.attackerPositions=new double[][]{')
					t=sim.topology.attackerPositions;
					for i=1:size(t,1)
						disp(['{',mat2str(t(i,1)),',',mat2str(t(i,2)),',',mat2str(t(i,3)),'},']);
					end
					disp('};')

					showTopology(sim.topology.nodePositions, ...
						sim.topologyParameters.communicationRange, ...
						sim.topology.sourceNodes+1, true);
				end

				%Initialize and run simulation
				sim.topology.createModel();
				sim.run();

				%Record results
				s_nodes(k)=sim.topology.nodeCount;
				sourceNodes=sim.topology.getSourceNodes();
				s_src_deg(k)=sim.topology.countSourceNeighbors(sourceNodes(1));
				s_genrtd(k)=sim.results.generatedCount;
				s_sent(k)=sim.results.sentCount; 
				s_recved(k)=sim.results.receivedCount;
				s_forwd(k)=sim.results.forwardedCount;
				s_pktcost(k)=(sim.results.forwardedCount/sim.results.receivedCount);
				s_e2edly(k)=sim.results.averageE2EDelay;
				s_fintim(k)=sim.results.finishTime;
				s_maxpktcost(k)=double(max(sim.results.numberOfTransmissions));
			end
			
			nodes(r)=round(100*mean(s_nodes))/100;
			src_deg(r)=round(100*mean(s_src_deg))/100;
			genrtd(r)=round(100*mean(s_genrtd))/100;
			sent(r)=round(100*mean(s_sent))/100; 
			recved(r)=round(100*mean(s_recved))/100;
			forwd(r)=round(100*mean(s_forwd))/100;
			pktcost(r)=round(100*mean(s_pktcost))/100;
			e2edly(r)=round(100*mean(s_e2edly(~isnan(s_e2edly)&(s_e2edly>0))))/100;
			maxpktcost(r)=round(100*mean(s_maxpktcost)/mean(s_recved))/100;
			fintim(r)=round(100*mean(s_fintim))/100;
			
		end
	end
	t= table(D, attacks, recved, forwd, pktcost,maxpktcost, e2edly,fintim, ...
	'VariableNames', {'D', 'attacks', 'recved', 'forwd', 'pktcost', 'maxpktcost', 'e2edly', 'fintim'});
	disp(t)
	results=[t.Properties.VariableNames;table2cell(t)];
end

function createTopology(sim, pos) 

    sim.topology.nodeCount=size(pos,1)-1; 
    sim.topology.setNodePositions(pos(:));
    sim.topology.seedRandom(sim.topologyParameters.topologySeed);
    sim.topology.sinkPos = pos(1,:);

    %sim.topology.chooseFurthestSource()
    for i=1:sim.topologyParameters.sourcesCount
        chooseDeepestSource(sim);		
%        chooseFurthestSource(sim)
    end
	sourceNodes=sim.topology.sourceNodes;
    for i=1:sim.topologyParameters.attackersCount 
        sim.topology.placeAttacker(sourceNodes(i));
    end
%    sim.topology.attackerPos = [271,68,84];

end


function chooseDeepestSource(sim) 
	pos=reshape(sim.topology.getNodePositions,[],3);
    depths=pos(:,3);
    depths(sim.topology.getSourceNodes()+1)=0;
    [~,i] = max(depths);
    sim.topology.addSourceNode(i-1);
end


function chooseFurthestSource(sim)
	pos=reshape(sim.topology.getNodePositions,[],3);
    distances=sqrt(sum(((ones(length(pos),1) * ...
        sim.topology.sinkPos') - pos).^2,2));
    distances(sim.topology.getSourceNodes()+1)=0;
    [~,i] = max(distances);
    sim.topology.addSourceNode(i-1);
end

