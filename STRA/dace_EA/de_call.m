function [bestvalue,bestind]=de_call(lb,ub,par)



nVar=size(lb,2);            % Number of Decision Variables设计变量维度
n=nVar;

VarSize=[1 nVar];   % Decision Variables Matrix Size 

lu=[lb;ub];

%% DE Parameters
if n<10
    MaxIt=100;      % Maximum Number of Iterations最大迭代次数
    nPop=50;        % Population Size种群大小
else
    MaxIt=80;      % Maximum Number of Iterations最大迭代次数
    nPop=100;        % Population Size种群大小
end


beta_min=0.2;   % Lower Bound of Scaling Factor  参数F的最小值
beta_max=0.8;   % Upper Bound of Scaling Factor  参数F的最大值
pCR_min=0.1;    % Lower Bound of Crossover Probability   参数CR的最小值 
pCR_max=0.5;    % Lower Bound of Crossover Probability   参数CR的最小值



%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop

    pop(i).Position=rand()*(lu(2,:)-lu(1,:))+lu(1,:);  %初始化种群个体的位置向量，这里采用的的是在设计变量的上下界范围内均匀分布的思想
    
    [pop(i).Cost,~]=objfunc(pop(i).Position, par);  %计算每个初始化的种群个体的适应度值
    
    if pop(i).Cost<BestSol.Cost
        BestSol.Cost=pop(i).Cost; %得到初始化种群中的最佳的位置向量。
        BestSol.ind=pop(i).Position;
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop

%设置DE的类型
detype=0;%detype=2时有问题，需要改进


FES=nPop;
for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];%这里是对第i个个体进行变异，因此，要先去掉第i个个体，然后随机。
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation  最关键的变异操作
        %beta=unifrnd(beta_min,beta_max);
        %这里beta是对于每一维度上的值都不一样，随机产生的，并不是一次迭代中所有个体都是一个定值，而是一个迭代中所有个体的设计变量的对应维度上的beta是一个定值，不同维度不一样，且不同迭代次数也不一样
        %因此，这里beta是一个1*Varsize的矩阵
        %这里设置参数
        switch detype 
            case 0
                %经典最简单的de,  F和CR均不变
                beta=0.8;
                pCR=0.8;        % Crossover Probability  交叉概率CR的值
            case 1
                %F动态变化，对每个维度上均不同；但是CR是固定的
                beta=unifrnd(beta_min,beta_max,VarSize);  %这里的缩放因子每次都取0.2-0.8范围内的一个随机数，当然也可以将其设置为一个定值。beta就是F
                pCR=0.2;        % Crossover Probability  交叉概率CR的值
            case 2
                %F动态变化，对每个维度上均不同；CR动态变化，对每个维度上均不同
                beta=unifrnd(beta_min,beta_max,VarSize);  %这里的缩放因子每次都取0.2-0.8范围内的一个随机数，当然也可以将其设置为一个定值。beta就是F
                pCR=unifrnd(pCR_min,pCR_max,VarSize);        % Crossover Probability  交叉概率CR的值
        end
        
        
%         switch mutation
%             case 1
%                 y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);%这里采用的就是简单的DE/rand/1的变异方法，也可以有其他的变异方法。要注意并不是最简单的普通的DE,其缩放因子beta在每一个维度上是不一样的值。
%             case 2
%                 y=BestSol.Position+
%         end
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);%这里采用的就是简单的DE/rand/1的变异方法，也可以有其他的变异方法。要注意并不是最简单的普通的DE,其缩放因子beta在每一个维度上是不一样的值。
        
        VarMin=repmat(lb,size(y,1),1);
        VarMax=repmat(ub,size(y,1),1);
        
        y = max(y, VarMin);
		y = min(y, VarMax);   %注意这个y是对应的变异向量，并不是对应于某个位置的目标函数的值。并且这里考虑到了如果变异向量超出去边界取边界！
        
        
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        [NewSol.Cost,~]=objfunc(NewSol.Position,par);
        FES=FES+1;
        

        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol.Cost=pop(i).Cost;
               BestSol.ind=pop(i).Position;
            end
        end
        
        
        
    end
    
    % Update Best Cost更新最佳适应度值
    BestCost(it)=BestSol.Cost;
    bestvalue=BestSol.Cost;
    bestind=BestSol.ind;
    
    % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end

%% Show Results


end
