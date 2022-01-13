clear
clc
close all

javarmpath('.\bin');
javaaddpath('.\bin');
import uansim.*

% load("pos_data.mat", "pos_data");
% pos_data=pos_data(1,:);
% pos_data=pos_data(:,1:1000);
% good_pos=false(size(pos_data));

load("good_pos_data.mat", "good_pos_data");
pos_data=good_pos_data(2:4,1:100);
dispTopology=false;
debugSimulation=false;
dispDetailedResults=false;
dispSummaryResults=true;
Networks=1:100;

algorithms={@updateForwardingProbability1,@updateForwardingProbability2,@updateForwardingProbabilityPolynomialFitting};
updateForwardingProbability=algorithms{2};
global polynomialDegree
polynomialDegree=4;
global eta1
global eta2
K=0.8;
eta1=K^2+0.1;
eta2=eta1;

sum_net=[];
sum_d=[];
sum_ps=[];
sum_rs=[];
sum_pc=[];
sum_rc=[];
sum_rptd=[];
sum_maxr=[];
sum_pktcost=[];
sum_pktcost=[];
sum_maxpktcost=[];

for net= Networks   %Network #
    pos_data=good_pos_data(2:4,net);   %Network #4, #9 polyfit fails when degree>3
                                    %Network #11, #12 attack failed

    %   network_density number_of_nodes min_number_of_neighbors
    experiments = [
    ...    4	28	3 1
        5	47	4 0.9
        8	75	7 0.5
        15	140	9 0.3
    ...    6	57	5 1
    ...    7	66	6 1
    ...    9	85	8 1
    ];

    %Result columns
    D=[]; p=[]; attacks=[]; nodes=[]; src_deg=[]; genrtd=[]; sent=[];
    recved=[]; forwd=[]; pktcost=[]; e2edly=[]; fintim=[]; maxpktcost=[];
    masterResults= table(D, p, attacks, recved, forwd, pktcost, maxpktcost, e2edly, fintim);

    %Max # of attackers
    A=1;

    % after each iteration, the forwarding probability is revised using the
    % algorithm defined by "updateForwardingProbability"

    iterations=30;

    for pi=1
        % for each density
        for ii = 1: size(pos_data,1)

            %for each initial forwarding probability
            for pp = 0 %:0.1:1.0 % pi*experiments(ii,4); %:0.1:1.0

                %for each number of attackers
                for a  = 1:A  
                    %history of forwarding probability
                    past=[];

                    %number of valid topologies
                    n=0;
                    %for each topology
                    for k= 1:size(pos_data,2)
                        %filter out topologies with invalid nodes
                        if any(all(pos_data{ii,k}==[0,0,0],2)) 
                            warning('Found topology with invalid node locations')
                        else
                            n=n+1;
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
                            sim.routingParameters.K = K;
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

                            if dispTopology && pi==1 && a==1 && pp==0 
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
                                set(gcf, 'NumberTitle', 'off', ...
                                    'Name', sprintf('Network #%d', net));
                                exportgraphics(gcf,['net',mat2str(experiments(ii,1)),'.png'],'Resolution',600)
                            end
                            d=experiments(mod(ii-1, length(experiments))+1,1);

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

                            % for each forwarding probability update period
                            for it=1:iterations
                                %Initialize and run simulation
                                sim.topology.createModel();
                                sim.run()

                                %Record results
                                D=[D; d];
                                attacks=[attacks;a];
                                s_nodes=[s_nodes;sim.topology.nodeCount];
                                s_p=[s_p;DPR.FORWARDING_PROBABILITY(sim.topology.sourceNodes+1)'];
                                s_src_deg=[s_src_deg;sim.topology.countSourceNeighbors(sim.topology.sourceNodes(1))];
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
                                % Update forwarding probability
                                past = updateForwardingProbability(sim, past,d);
                            end
                        end
    %                     nodes=[mean(s_nodes);sim.topology.nodeCount];
    %                     p=[p;round(100*mean(s_p),2)];                
    %                     src_deg=[src_deg;round(mean(s_src_deg),2)];
    %                     genrtd=[genrtd;round(mean(s_genrtd),2)];
    %                     sent=[sent;round(mean(s_sent),2)]; 
    %                     recved=[recved;round(mean(s_recved),2)];
    %                     forwd=[forwd;round(mean(s_forwd),2)];
    %                     pktcost=[pktcost;round(sum(s_forwd)/sum(s_recved),2)];                
    %                     e2edly=[e2edly;round(mean(s_e2edly(~isnan(s_e2edly)&(s_e2edly>0))),2)];
    %                     fintim=[fintim;round(mean(s_fintim(~isnan(s_fintim))),2)];
    %                     maxpktcost=[maxpktcost;round(mean(s_maxpktcost)/mean(s_recved),2)];                
                        nodes=[s_nodes;sim.topology.nodeCount];
                        p=[p;s_p];                
                        src_deg=[src_deg;s_src_deg];
                        genrtd=[genrtd;s_genrtd];
                        sent=[sent;s_sent]; 
                        recved=[recved;s_recved];
                        forwd=[forwd;s_forwd];
                        pktcost=[pktcost;s_forwd./s_recved];                
                        e2edly=[e2edly;s_e2edly];
                        fintim=[fintim;s_fintim];
                        maxpktcost=[maxpktcost;s_maxpktcost./s_recved];                
                    end
    %                disp(['Number of valid topologies = ' mat2str(n)])
                    % K=sim.routingParameters.K * ones(iterations,1);
                    results= table(D, p, attacks, recved, forwd, pktcost, maxpktcost, e2edly, fintim);

                    if dispDetailedResults
                        disp(results)
                    end

                    if dispSummaryResults
                        disp(['Net #' mat2str(net), ...
                            '   , density=', mat2str(d) ...
                            '   , earliest convergence=' mat2str(find(recved>=100*K-2,1))...
                            '   , number of successes=' mat2str(sum(recved>=100*K-2))...
                            '   , converged at p=' mat2str(p(end)) ...
                            '   and r=' mat2str(recved(end)) ...
                            '   repeated=' mat2str(sum(p==p(end)))...
                            '   max r=' mat2str(max(recved))])
                    end
                    sum_net=[sum_net;net];
                    sum_d=[sum_d;d];
                    sum_ps=[sum_ps;p'];
                    sum_rs=[sum_rs;recved'];
                    sum_pc=[sum_pc;p(end)];
                    sum_rc=[sum_rc;recved(end)];
                    sum_rptd=[sum_rptd;sum(p==p(end))];
                    sum_maxr=[sum_maxr;max(recved)];
                    sum_pktcost=[sum_pktcost;pktcost];
                    sum_maxpktcost=[sum_maxpktcost;maxpktcost];
                    
                    masterResults=[masterResults;results];  
                    D=[]; p=[]; attacks=[]; nodes=[]; src_deg=[]; genrtd=[]; sent=[];
                    recved=[]; forwd=[]; pktcost=[]; e2edly=[]; fintim=[]; maxpktcost=[];
                end
            end
        end
    end
    % save("good_pos.mat", "good_pos");
end

results=table(sum_net,sum_d,sum_ps,sum_rs,sum_pc,sum_rc,sum_rptd,sum_maxr)


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


function createTopology1(sim, nodeCount) 

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

function past=updateForwardingProbability1(sim, past, d)
    if isempty(past)
        past=zeros(sim.topologyParameters.sourcesCount,1);
    end
    p = sim.routingParameters.forwardingProbability;
    r = sim.results.deliveryRatio;
    K=sim.routingParameters.K;
    ts=past;
    for i=1:length(sim.topology.sourceNodes)
        j=sim.topology.sourceNodes(i)+1;
        rs=r(i);
        if abs(rs-K)>0.02
            ts(i)=K*(ts(i)-rs+1);
            if ts(i)>1
                ts(i)=1;
            end
            p(j)=min(1.0,round(100*(1-nthroot(1-ts(i),d)))/100);
        end
    end
    sim.routingParameters.forwardingProbability=p;
    past=ts;
end

function past=updateForwardingProbability2(sim, past, d)
    global eta1
    global eta2

    if isempty(past)
        past=zeros(sim.topologyParameters.sourcesCount,1);
    end
    p = sim.routingParameters.forwardingProbability;
    r = sim.results.deliveryRatio;
   past = [r, past];
    K=sim.routingParameters.K;
    for i=1:length(sim.topology.sourceNodes)
        j=sim.topology.sourceNodes(i)+1;
%         ts=past(i);
        ts=p(j);
        rs=r(i);
        rs=max(r(i),mean(past(i,1:min([30,min(size(past,2))]))));
        if K-rs > 0.02 
            %Increase threat level
            ts2 = (1 - eta1*(K - rs))*ts + eta1*(K - rs);
        elseif K-rs < -0.02 && ts>0 
            %Decrease threat level
            ts2 = max([0, ts + eta2 * (K - rs)]);
        elseif ts<0
            ts2=0;
        else
            ts2=ts;
        end
%         p(j)=1-nthroot(1-ts2,d);
        p(j) = min(1.0,round(100*ts2)/100);
%         past(i)=ts2;
    end
    sim.routingParameters.forwardingProbability=p;
end

function [past]=updateForwardingProbabilityPolynomialFitting(sim, past, d)
    global polynomialDegree
    if isempty(past)
        past_p={};
        past_r={};
        past_p{1:sim.topologyParameters.sourcesCount}=[]; 
        past_r{1:sim.topologyParameters.sourcesCount}=[];
    else
        %past_p for sources only
        past_p = past{1};
        %past_r for sources only
        past_r = past{2};    
    end

%
    %Estimated number of unqualified forwarders
    D=100;
    
    %Resolution of p = 1/R
    p_R=100;
    
    %Resolution of r
    r_R=100;
    
    %%target delivery ratio
    K=sim.routingParameters.K;

    p = sim.routingParameters.forwardingProbability;
    r = round(sim.results.deliveryRatio*r_R)/r_R;
    
    for i=1:length(sim.topology.sourceNodes)
        j=sim.topology.sourceNodes(i)+1;

        %remove points that are too close to the new point
        p_p=past_p{i};
        p_r=past_r{i};
        p_r=p_r(abs(p_p-p(j))>0.01);        %Clean the series from the points that are too close to the most recent point
        p_p=p_p(abs(p_p-p(j))>0.01);
        
        %add the new points
        p_p = [p(j), p_p];   
        p_r = [r(i), p_r];
%         p_q = p_r./(p_p+1/D);  %estimate the network efficiency
        p_q = p_r;  %estimate the network efficiency
        past_p{i}=p_p;
        past_r{i}=p_r;

        n = min(length(p_p)-1,polynomialDegree);
    
        if n==0 % && r(i)<K 
            p(j)=0.5;
        elseif n==1 % && r(i)<K 
            p(j)=0.75;
        elseif n==2 % && r(i)<K
            p(j)=0.25;
        else
            poly = polyfit(p_p, p_r, n);
            if false
                figure
                plot((0:0.05:1),polyval(poly,(0:0.05:1)))
                hold on 
                plot(p_p,p_r,'x')
                title(['n=',mat2str(n)])
            end
            %%solve poly-K==0
            xs=roots([poly(1:n),poly(n+1)-K]);
            %%solution must be in the interval [0,1]
            xs=xs(real(xs)==xs & xs>=0 & xs<=1);

            if ~isempty(xs)     %if solutions exists
                %p(j)=min(xs);   %take the smallest one
                p(j)=mean([p(j),min(xs)]);   %move towards the smallest one
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
        p(j)=round(p_R*p(j))/p_R;
    end
    sim.routingParameters.forwardingProbability=p;
    past{1}=past_p;
    past{2}=past_r;
end

