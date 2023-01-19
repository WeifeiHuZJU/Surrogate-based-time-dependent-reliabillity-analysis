function [cell]=voronoi(sam_domain_x,sam_x,cellposition)


n=size(sam_domain_x,1);
m=size(sam_x,1);
cell=[];


for i=1:n
    distance=zeros(m,1);
    %计算第i个候选点到所有已有点的距离（二阶范数）
    for j=1:m
        distance(j)=norm(sam_domain_x(i,:)-sam_x(j,:));
    end
    [~,pos]=min(distance);%找到距离最小的点
    
    %判断第i个点是否属于最敏感的cell
    if(ismember(pos,cellposition))
        cell=[cell;sam_domain_x(i,:)];
    end
end

end