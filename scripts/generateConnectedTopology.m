function [pos]=generateConnectedTopology(sim)
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
            maxTrials=10000;
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

