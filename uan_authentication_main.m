clear
clc
close all

javaaddpath('.\bin');
import uansim.*

debug=true;
dispTopology=false;


netload=[]; msglen=[]; attacks=[];nodes=[]; src_deg=[];genrtd=[];sent=[];recved=[];forwd=[];pktcost=[];
past=[];

overhead=160;
for genrate=[10,50,100]
    for msglength = [100,1000]
        for atkCount=0:7 % for scenario #2
            sim=Simulator();
            sim.debug = debug;

            %%%%%%%%%% PARAMETERS %%%%%%%%%%
            sim.topologyParameters.topologySeed = 17943;
            sim.topologyParameters.networkDensity = 3;
            sim.topologyParameters.deploymentSideLength=500;
            sim.topologyParameters.deploymentDepth=250;
            sim.topologyParameters.sourcesCount = 7;
            sim.topologyParameters.attackersCount =atkCount; %(change in scenario 2,3)

            sim.topologyParameters.communicationRange=150;
            sim.topologyParameters.interferenceRange=2*150;
            nodeCount=30;

            sim.trafficParameters.messageCount=100;      %(per source)
            sim.trafficParameters.generationBitRate=genrate; %(to control load)
            sim.trafficParameters.messageBitLength=msglength;
            sim.trafficParameters.channelBitRate=10000;
            sim.trafficParameters.framingBitOverhead=overhead; %(change in scenario 3,4)

            sim.routingParameters.K = 0.0;
            sim.routingParameters.forwardingProbability = 0.0*ones(256,1);
            sim.routingParameters.forwardingThreshold = 0;
            sim.routingParameters.maximumForwardingDelay=25;

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
                    sim.topology.sourceNodes+1);
            end

            sim.topology.createNodes();
            sim.run()
            netload=[netload; 100*sim.trafficParameters.generationBitRate*sim.topologyParameters.sourcesCount / (sim.trafficParameters.channelBitRate*0.09)];
            msglen=[msglen; sim.trafficParameters.messageBitLength];
            attacks=[attacks;sim.topologyParameters.attackersCount];
            nodes=[nodes;sim.topology.nodeCount];
            src_deg=[src_deg;sim.topology.countSourceNeighbors(sim.topology.sourceNodes(1))];
            genrtd=[genrtd;sim.results.generatedCount];
            sent=[sent;sim.results.sentCount]; 
            recved=[recved;sim.results.receivedCount];
            forwd=[forwd;sim.results.forwardedCount];
            pktcost=[pktcost;round(sim.results.forwardedCount/sim.results.receivedCount,2)];

            %    sim.updateForwardingProbability();
            %    past = updateForwardingProbability(sim, past);
        end
    end
end
results= table(netload,msglen, attacks, nodes, genrtd, recved, round(100*recved./genrtd,1), pktcost);
disp(results)

function createTopology(sim, nodeCount) 

    sim.topology.nodeCount=nodeCount;
    
    generateConnectedTopology(sim);
%    generateRegularTopology(sim);

    %sim.topology.chooseFurthestSource()
    for i=1:sim.topologyParameters.sourcesCount
        chooseDeepestSource(sim);		
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
    for i=2:n+1
        bad=true;
        while bad
            bad = false;
            p = top.randVector();
            newNodeDegree=nodeDegree;
            for j=1:i-1
                v = pos(j,:);
                if top.inCommunicationRange(p,v)
                    if newNodeDegree(j) > params.networkDensity || ...
                        newNodeDegree(i) > params.networkDensity
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
        end

        nodeDegree = newNodeDegree;
        pos(i,:)=p;
    end
    
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
        else
            ts2 = ts + 0.5 * (K - rs);
        end
        p(j) = ts2;
    end
    sim.routingParameters.forwardingProbability=p;
end
