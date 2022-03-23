%% Dynamic Time Warping (DTW)

function [dtw_dist, dtw_matrix]=DTW(x1,x2)
x1_n=size(x1,1);
x2_n=size(x2,1);

% plot(x1,'b');
% hold on
% plot(x2,'r');

distance=(repmat(x1(:),1,x1_n)-repmat(x2(:)',x2_n,1)).^2;

cost=zeros(x1_n+1,x2_n+1);
cost(:,1)=inf;
cost(1,:)=inf;
cost(1,1)=0;

for i=1:x1_n
    for j=1:x2_n
        cost(i+1,j+1)=distance(i,j)+min([cost(i,j),cost(i,j+1),cost(i+1,j)]);
    end
end
dtw_matrix=cost(2:end,2:end);
dtw_dist=cost(x1_n+1,x2_n+1);
end