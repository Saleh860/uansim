clear
clc
close all

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
deployments = 10;
pos_data={};

% for each experiment
for ii = 1: size(experiments,1) 
        %for each number of attackers
        jj=1;
        while jj<=deployments
            sim=Simulator();

            %%%%%%%%%% PARAMETERS %%%%%%%%%%
            % Topology parameters
            sim.topologyParameters.topologySeed = randi(30000);
            sim.topologyParameters.networkDensity = experiments(mod(ii-1, size(experiments,1))+1,1);
            sim.topologyParameters.deploymentSideLength=500;
            sim.topologyParameters.deploymentDepth=250;
            sim.topologyParameters.sourcesCount = 1;
            sim.topology.nodeCount = experiments(mod(ii-1, length(experiments))+1,2);

            dispTopology=false;
            %%%%%%%%%% Create nodes %%%%%%%%%%%
            %sim.topology.prepare();
            
            pos=generateConnectedTopology(sim);
            disp(length(pos))
            if(length(pos)>sim.topology.nodeCount)
                pos_data{ii,jj}=pos;
                jj=jj+1;
            end
            
        end
end

save("pos_data.mat", "pos_data");
