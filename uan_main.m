clear
clc
%close all

javaaddpath('.\bin');
import uansim.*

debug=false;
dispTopology=false;
sourcesCount=2;

experiments = [
    4	28	3
...    5	47	4
...    6	57	5
...    7	66	6
...    8	75	7
...    9	85	8
...    15	140	9
];

it=[]; p=[]; attacks=[];nodes=[]; src_deg=[];genrtd=[];sent=[];recved=[];
Ks=[]; pdr=[]; forwd=[];pktcost=[];
past_p={};
past_r={};
for i=1:sourcesCount
    past_p{i}=[]; past_r{i}=[];
end

for ii = 1: size(experiments,1) 
    for a  = 0:sourcesCount %1  
        for K = 0:0.1:1.0

            sim=Simulator();


            %%%%%%%%%% PARAMETERS %%%%%%%%%%
            sim.topologyParameters.topologySeed = 28629;
            sim.topologyParameters.networkDensity = experiments(mod(ii-1, size(experiments,1))+1,1);
            sim.topologyParameters.deploymentSideLength=500;
            sim.topologyParameters.deploymentDepth=250;
            sim.topologyParameters.sourcesCount = sourcesCount;
            sim.topologyParameters.attackersCount = a;
            nodeCount = experiments(mod(ii-1, length(experiments))+1,2);

            sim.trafficParameters.messageCount=100;
            sim.trafficParameters.generationRate=10;
            sim.trafficParameters.messageLength=50;
            sim.trafficParameters.channelBitRate=100000/8;

            sim.routingParameters.K = K;
            sim.routingParameters.forwardingProbability = 0.0*ones(256,1);
            sim.routingParameters.forwardingThreshold = 0;
            sim.routingParameters.maximumForwardingDelay=1;

            iterations=100;
            decisionSeeds=[36411,14283,29685,18527,39782,3677,49637,38802,369,44448,43207,34880,44250,23283,18578,4152,49861,43492,4,33999,32745,20591,23819,19414,22603,1338,21265,5291,2156,27807,30565,36942,12846,36313,35663,25189,44639,38970,24744,38076,33353,48906,5854,43467,28393,6008,48517,41999,3116,13904,20482,31341,45956,20758,35576,48544,21774,31770,45078,43580,7136,10030,31047,38752,12341,18473,17055,22820,8542,37281,4899,44514,4105,39077,42434,22544,30267,17724,37955,29473,30579,42159,22850,44140,24679,31143,28109,36084,27339,39223,22318,45596,4330,49947,36380,27921,7142,41882,36263,20869];
            trafficSeeds=decisionSeeds+1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            sim.debug = debug;
            %sim.topology.prepare();
            createTopology(sim, nodeCount);
            for i=1:iterations
                sim.routingParameters.decisionSeed=decisionSeeds(mod(i-1,length(decisionSeeds))+1);
                sim.trafficParameters.trafficSeed=trafficSeeds(mod(i-1,length(trafficSeeds))+1);
                sim.topology.createNodes();
                sim.run()

                p=[p;round(100*DPR.FORWARDING_PROBABILITY(sim.topology.sourceNodes+1))'];
                attacks=[attacks;sim.topologyParameters.attackersCount];
                nodes=[nodes;sim.topology.nodeCount];
                src_deg=[src_deg;sim.topology.countSourceNeighbors(sim.topology.sourceNodes(1))];
                genrtd=[genrtd;sim.results.generatedCount];
                sent=[sent;sim.results.sentCount]; 
                recved=[recved;sim.results.receivedCount];
                pdr=[pdr;round(100*sim.results.deliveryRatio')];
                forwd=[forwd;sim.results.forwardedCount];
                pktcost=[pktcost;round(sim.results.forwardedCount/sim.results.receivedCount,2)];

            %    sim.updateForwardingProbability();
                [past_p, past_r] = updateForwardingProbability(sim, past_p, past_r);
            end
            it=[it;(1:iterations)'];
            Ks=[Ks;round(100*sim.routingParameters.K) * ones(iterations,1)];            
        end    
    end

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

    end

    if dispTopology
        showTopology(sim.topology.nodePositions, ...
            sim.topologyParameters.communicationRange, ...
            sim.topology.sourceNodes+1);
    end

end
results= table(it, p, Ks, attacks, nodes, src_deg, genrtd, sent, recved, pdr, forwd, pktcost);
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
    for i=2:n+1
        bad=true;
        while bad
            bad = false;
            p = top.randVector();
            newNodeDegree=nodeDegree;
            for j=1:i-1
                v = pos(j,:);
                if top.inCommunicationRange(p,v)
                    if newNodeDegree(j) >= params.networkDensity-1 || ...
                        newNodeDegree(i) >= params.networkDensity-1
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

function [past_p,past_r]=updateForwardingProbability(sim, past_p, past_r)
%past_p for sources only
%past_r for sources only
%
    %Estimated number of unqualified forwarders
    D=100;
    
    %Resolution of p = 1/R
    R1=100;
    
    %Resolution of r
    R2=100;
    
    %%target delivery ratio
    K=sim.routingParameters.K;

    p = sim.routingParameters.forwardingProbability;
    r = round(sim.results.deliveryRatio*R2)/R2;
    
    for i=1:length(sim.topology.sourceNodes)
        j=sim.topology.sourceNodes(i)+1;

        %remove points that are too close to the new point
        p_p=past_p{i};
        p_r=past_r{i};
        p_r=p_r(abs(p_p-p(j))>0.01);
        p_p=p_p(abs(p_p-p(j))>0.01);
        
        %add the new points
        p_p = [p(j), p_p];   
        p_r = [r(i), p_r];
        p_q = p_r./(p(j)+1/D);
        past_p{i}=p_p;
        past_r{i}=p_r;

        %% polynomial degree 2
        n = min(length(p_p)-1,2);
    
        if r(i)<K && n==0
            p(j)=0.5;
        else
            poly = polyfit(p_p, p_r, n);
            %%solve poly-K==0
            xs=roots([poly(1:n),poly(n+1)-K]);
            %%solution must be in the interval [0,1]
            xs=xs(real(xs)==xs & xs>=0 & xs<=1);

            if length(xs)>0     %if solutions exists
                p(j)=min(xs);   %take the smallest one
            else    %maximize delivery ratio
                poly = polyfit(p_p, p_q, n);
                poly1=polyder(poly);
                poly2=polyder(poly1);
                xs=roots(poly1)';
                xs=xs(real(xs)==xs);
                xs=[0,1,xs(polyval(poly2,xs)<0 & xs>=0 & xs<=1)];
                ys=polyval(poly,xs);
                p(j)=min(xs(ys==max(ys)));
            end
        end
        p(j)=round(R1*p(j))/R1;
    end
    sim.routingParameters.forwardingProbability=p;
end
