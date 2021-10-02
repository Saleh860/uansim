clear
clc
close all

% load("pos_data.mat", "pos_data");
% pos_data=pos_data(1,:);
% pos_data=pos_data(:,1:10);

load("good_pos_data.mat", "good_pos_data");
pos_data=good_pos_data(:,1:100);

javarmpath('.\bin');
javaaddpath('.\bin');
import uansim.*

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
            sim=RPRSimulator();

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
            sim.debug = false;
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
            sim.run()

            %Record results
            s_nodes(k)=sim.topology.nodeCount;
            s_src_deg(k)=sim.topology.countSourceNeighbors(sim.topology.sourceNodes(1));
            s_genrtd(k)=sim.results.generatedCount;
            s_sent(k)=sim.results.sentCount; 
            s_recved(k)=sim.results.receivedCount;
            s_forwd(k)=sim.results.forwardedCount;
            s_pktcost(k)=(sim.results.forwardedCount/sim.results.receivedCount);
            s_e2edly(k)=sim.results.averageE2EDelay;
            s_fintim(k)=sim.results.finishTime;
            s_maxpktcost(k)=double(max(sim.results.numberOfTransmissions));
        end
        
        nodes(r)=round(mean(s_nodes),2);
        src_deg(r)=round(mean(s_src_deg),2);
        genrtd(r)=round(mean(s_genrtd),2);
        sent(r)=round(mean(s_sent), 2); 
        recved(r)=round(mean(s_recved),2);
        forwd(r)=round(mean(s_forwd),2);
        pktcost(r)=round(mean(s_pktcost),2);
        e2edly(r)=round(mean(s_e2edly(~isnan(s_e2edly)&(s_e2edly>0))),2);
        maxpktcost(r)=round(mean(s_maxpktcost)/mean(s_recved),2);
        fintim(r)=round(mean(s_fintim),2);
        
    end
end
results= table(D, attacks, recved, forwd, pktcost,maxpktcost, e2edly,fintim);
disp(results)

function createTopology(sim, pos) 

    sim.topology.nodeCount=size(pos,1)-1; 
    sim.topology.nodePositions=pos;
    sim.topology.seedRandom(sim.topologyParameters.topologySeed);
    sim.topology.sinkPos = pos(1,:);

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

