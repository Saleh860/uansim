function showTopology(pos, range, sources, showNodeLabels)
% pos=[
%     0   0 0
%     0   100 20
%     0   200 30
%     100 200 40
%     200 0   10
%     ];
% range = 150;
%%
n=size(pos,1);
adjacency=range>=sqrt((pos(:,1)*ones(1,n) - ones(n,1)*pos(:,1)').^2 + (pos(:,2)*ones(1,n) - ones(n,1)*pos(:,2)').^2 + (pos(:,3)*ones(1,n) - ones(n,1)*pos(:,3)').^2);
%%
f=figure;
if nargin<4
    showNodeLabels=false;
end
if nargin==2
    plot3(pos(:,1), pos(:,2), -pos(:,3), '.', 'MarkerSize',15);
    hold on
else
    plot3(pos(sources,1), pos(sources,2), -pos(sources,3), 'r.', 'MarkerSize',30);
    hold on
    plot3(pos(1,1), pos(1,2), -pos(1,3), 'g.', 'MarkerSize',20);
    hold on
    ii=true(n,1);
    ii(1)=false;
    ii(sources)=false;
    plot3(pos(ii,1), pos(ii,2), -pos(ii,3), 'b.', 'MarkerSize',20);
    xlabel('(m)')
    ylabel('(m)')
    zlabel('(m)')
    if showNodeLabels
        for i=1:length(pos)
            text(pos(i,1),pos(i,2),-pos(i,3),num2str(i-1), 'fontsize',12); 
        end
    end
end
    
for i=1:n
    for j=1:n
        if adjacency(i,j)
            plot3(pos([i,j],1), pos([i,j],2), -pos([i,j],3), '--');
        end
    end
end
f.Position = [1 1 720 720];
set(gca,'FontSize',14)