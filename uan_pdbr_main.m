clear
clc
close all

load("pos_data.mat", "pos_data");
% pos_data=pos_data(1,:);

javarmpath('.\bin');
javaaddpath('.\bin');
import uansim.*

%   network_density number_of_nodes min_number_of_neighbors
experiments = [
...    4	28	3
    5	47	4
...    6	57	5
...    7	66	6
    8	75	7
...    9	85	8
    15	140	9
];

%Result columns
D=[]; p=[]; attacks=[];nodes=[]; src_deg=[];genrtd=[];sent=[];recved=[];forwd=[];pktcost=[];


%Max # of attackers
A=1;

%after each iteration, the forwarding probability is revised
iterations=10;

% for each experiment
for ii = 1: size(pos_data,1) 
    %for each initial forwarding probability
    for pp = 0:0.1:1.0
        %for each number of attackers
        for a  = 0:1  

            sim=Simulator();

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

            %Simulation parameters
            sim.debug = false;
            dispTopology=false;
            %%%%%%%%%% Create nodes %%%%%%%%%%%
            %sim.topology.prepare();
            createTopology(sim, nodeCount);
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

            %history of forwarding probability
            past=[];
            
            % for each forwarding probability update period
            for i=1:iterations
                %Initialize and run simulation
                sim.topology.createModel();
                sim.run()
                
                %Record results
                D=[D; experiments(mod(ii-1, length(experiments))+1,1)];
                p=[p;round(100*DPR.FORWARDING_PROBABILITY(sim.topology.sourceNodes+1))'];
                attacks=[attacks;sim.topologyParameters.attackersCount];
                nodes=[nodes;sim.topology.nodeCount];
                src_deg=[src_deg;sim.topology.countSourceNeighbors(sim.topology.sourceNodes(1))];
                genrtd=[genrtd;sim.results.generatedCount];
                sent=[sent;sim.results.sentCount]; 
                recved=[recved;sim.results.receivedCount];
                forwd=[forwd;sim.results.forwardedCount];
                pktcost=[pktcost;round(sim.results.forwardedCount/sim.results.receivedCount,2)];

                % Update forwarding probability
%                 past = updateForwardingProbability(sim, past);
            end
        end
    end 
end
K=sim.routingParameters.K * ones(iterations,1);
results= table(D, p, attacks, recved, forwd, pktcost);
disp(results)

function createTopology(sim, nodeCount) 

    sim.topology.nodeCount=nodeCount;
    
    generateConnectedTopology(sim);
%    generateRegularTopology(sim);

    %sim.topology.chooseFurthestSource()
    for i=1:sim.topologyParameters.sourcesCount
        chooseDeepestSource(sim);		
%        chooseFurthestSource(sim)
    end
    for i=1:sim.topologyParameters.attackersCount 
        sim.topology.placeAttacker(sim.topology.sourceNodes(i));
    end
%    sim.topology.attackerPos = [271,68,84];

end

function generateConnectedTopology(sim)
    top=sim.topology;
    params=sim.topologyParameters;
    n = sim.topology.nodeCount;
    top.seedRandom(params.topologySeed);
	% Good nodes locations, including the sink
    pos=zeros(n+1,3);
    top.sinkPos = params.deploymentSideLength/2*[1,1,0];
    pos(1,:) = top.sinkPos;
    nodeDegree = zeros(n+1,1);
    startNode = 2;
    bad = true;
    for currentNodeDegree = 3 : params.networkDensity-1
        for i=startNode:n+1
            bad=true;
            numTrials=0;
            maxTrials=1000;
            while bad && numTrials < maxTrials
                bad = false;
                p = top.randVector();
                newNodeDegree=nodeDegree;
                for j=1:i-1
                    v = pos(j,:);
                    if top.inCommunicationRange(p,v)
                        if newNodeDegree(j) >= currentNodeDegree || ...
                            newNodeDegree(i) >= currentNodeDegree
                            bad=true;
                            break;
                        else
                            newNodeDegree(i) = newNodeDegree(i)+1;
                            newNodeDegree(j) = newNodeDegree(j)+1;
                        end
                    end
                end
                if newNodeDegree(i)<1
                    bad=true;
                end
                if ~bad
                    bad=true;
                    for j=1:i-1
                        v = pos(j,:);
                        if top.inCommunicationRange(v,p) && v(3)<p(3)
                            bad=false;
                            break;
                        end
                    end
                end
                numTrials=numTrials+1;
            end
            if ~bad
                nodeDegree = newNodeDegree;
                pos(i,:)=p;
                startNode = i+1;
            else
                startNode = i;
                break;
            end
        end
    end
    if bad
        disp(['Actual number of nodes = ' mat2str(startNode-2)]);
        sim.topology.nodeCount = startNode-2;
        n = sim.topology.nodeCount;
    end
%    disp(['Maximum node degree = ' , num2str(max(nodeDegree))])
    for i=1:n+1
        for j=1:n+1
            adjacency(i,j) = top.inCommunicationRange(...
                pos(i,:), pos(j,:));
        end
    end

    top.nodePositions=pos;
end

function chooseDeepestSource(sim) 
    depths=sim.topology.nodePositions(:,3);
    depths(sim.topology.sourceNodes+1)=0;
    [~,i] = max(depths);
    sim.topology.sourceNodes=[sim.topology.sourceNodes; i-1];
end

function chooseFurthestSource(sim)
    distances=sqrt(sum(((ones(length(sim.topology.nodePositions),1) * ...
        sim.topology.sinkPos') - sim.topology.nodePositions).^2,2));
    distances(sim.topology.sourceNodes+1)=0;
    [~,i] = max(distances);
    sim.topology.sourceNodes=[sim.topology.sourceNodes; i-1];
end


function past=updateForwardingProbability(sim, past)
    p = sim.routingParameters.forwardingProbability;
    r = sim.results.deliveryRatio;
    past = [r, past];
    K=sim.routingParameters.K;
    for i=1:length(sim.topology.sourceNodes)
        j=sim.topology.sourceNodes(i)+1;
        rs=max(r(i),mean(past(i,1:min([30,min(size(past,2))]))));
        ts=p(j);
        if rs<=K
            ts2 = (1 + rs - K)*ts + (K - rs);
        elseif ts>0
            ts2 = ts + 0.5 * (K - rs);
        else
            ts2=0;
        end
        p(j) = ts2;
    end
    sim.routingParameters.forwardingProbability=p;
end

