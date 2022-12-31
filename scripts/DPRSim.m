function [results,t]=DPRSim(basedir,datadir)
	
	dispTopology=false;
	debugSimulation=0;

	datadir=config(basedir,datadir);	
	load(fullfile(datadir, 'data','good_pos_data.mat'), "good_pos_data");
	pos_data=good_pos_data(2:4,1:100);

	%   network_density number_of_nodes min_number_of_neighbors
	experiments = [
	...    4	28	3 1
		5	47	4 0.9
	...    6	57	5 1
	...    7	66	6 1
		8	75	7 0.5
	...    9	85	8 1
		15	140	9 0.3
	];

	%Result columns
	D=[]; p={}; attacks=[]; nodes=[]; src_deg={}; genrtd={}; sent={};
	recved={}; forwd={}; pktcost={}; e2edly={}; fintim={}; maxpktcost={};

	%Max # of attackers
	A=1;

	%after each iteration, the forwarding probability is revised
	iterations=1;

	for pi=1
		% for each experiment
		for ii = 1: size(pos_data,1) 
			%for each initial forwarding probability
			for pp = 0:0.1:1.0 % pi*experiments(ii,4); %:0.1:1.0
				%for each number of attackers
				for a  = 0:A  
					%history of forwarding probability
					past=[];

					% for each forwarding probability update period
					for it=1:iterations

						D=[D; experiments(mod(ii-1, length(experiments))+1,1)];
						attacks=[attacks;a];

						s_p=[];
						s_nodes=[];
						s_src_deg=[];
						s_genrtd=[];
						s_sent=[]; 
						s_recved=[];
						s_forwd=[];
						s_pktcost=[];
						s_e2edly=[];
						s_maxpktcost=[];
						s_fintim=[];

						%number of valid topologies
						n=0;
						for k= 1:size(pos_data,2)
							%filter out topologies with invalid nodes
							if any(all(pos_data{ii,k}==[0,0,0],2)) 
								warning('Found topology with invalid node locations')
							else
								n=n+1;
								sim=createJavaObject('uansim.Simulator');

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
								sim.trafficParameters.generationBitRate=10*8;
								sim.trafficParameters.messageBitLength=50*8;
								sim.trafficParameters.channelBitRate=100000/8;

								% Routing parameters
								sim.routingParameters.K = .0;
								sim.routingParameters.forwardingProbability = pp*ones(256,1);
								sim.routingParameters.forwardingThreshold = 0;
								sim.routingParameters.maximumForwardingDelay=1;
								sim.routingParameters.decisionSeed = 64641;
								sim.routingParameters.maximumTTL=20;

								%Simulation parameters
								sim.debug = debugSimulation;
								%%%%%%%%%% Create nodes %%%%%%%%%%%
								%sim.topology.prepare();
				%                createTopology(sim, nodeCount);
								%transfer topology to simulator
								createTopology(sim, pos_data{ii,k});

								if dispTopology && pi==1 && a==0 && pp==0 && it==1
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
										sim.topology.sourceNodes+1, false);
									exportgraphics(gcf,['net',mat2str(experiments(ii,1)),'.png'],'Resolution',600)
								end
								%Initialize and run simulation
								sim.topology.createModel();
								sim.run();

								%Record results
								s_nodes=[s_nodes;sim.topology.nodeCount];
								fp=sim.routingParameters.forwardingProbability;
								sourceNodes=sim.topology.getSourceNodes();
								s_p=[s_p;fp(sourceNodes+1)'];
								s_src_deg=[s_src_deg;sim.topology.countSourceNeighbors(sourceNodes(1))];
								s_genrtd=[s_genrtd;sim.results.generatedCount];
								s_sent=[s_sent;sim.results.sentCount]; 
								s_recved=[s_recved;sim.results.receivedCount];
								s_forwd=[s_forwd;sim.results.forwardedCount];
								s_e2edly=[s_e2edly;sim.results.averageE2EDelay];
								s_fintim=[s_fintim;sim.results.finishTime];
								s_maxpktcost(k)=double(max(sim.results.numberOfTransmissions));
								if pi==0 && a==0 && sim.results.receivedCount==sim.results.generatedCount
									good_pos(ii,k)=true;
								end
							end
						end
						nodes=[mean(s_nodes);sim.topology.nodeCount];
						p=[p;sprintf('%6.2f',100*mean(s_p))];                
						src_deg=[src_deg;sprintf('%6.2f',mean(s_src_deg))];
						genrtd=[genrtd;sprintf('%6.2f',mean(s_genrtd))];
						sent=[sent;sprintf('%6.2f',mean(s_sent))]; 
						recved=[recved;sprintf('%6.2f',mean(s_recved))];
						forwd=[forwd;sprintf('%6.2f',mean(s_forwd))];
						pktcost=[pktcost;sprintf('%6.2f',sum(s_forwd)/sum(s_recved))];
						e2edly=[e2edly;sprintf('%6.2f',mean(s_e2edly(~isnan(s_e2edly)&(s_e2edly>0))))];
						fintim=[fintim;sprintf('%6.2f',mean(s_fintim(~isnan(s_fintim))))];
						maxpktcost=[maxpktcost;sprintf('%6.2f',mean(s_maxpktcost)/mean(s_recved))];
					end
				end

			end
		end
	end
	% disp(['Number of valid topologies = ' mat2str(n)])
	K=sim.routingParameters.K * ones(iterations,1);
	t = table(		D, 	p, 	attacks, recved, forwd, 		pktcost, 	maxpktcost, 	e2edly, 	fintim, ...
	'VariableNames',{	'D','p','atck#', 'r'   , 'total_xmits', 'xmits/pkt', 'max_xmits/pkt', 'E2E_delay', 'finish_time'});
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

	sourceNodes=sim.topology.getSourceNodes();
    for i=1:min([sim.topologyParameters.attackersCount , length(sourceNodes)])
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



