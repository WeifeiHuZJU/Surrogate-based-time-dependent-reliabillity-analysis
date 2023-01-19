function [bestvalue,bestind]=de_call(lb,ub,par)



nVar=size(lb,2);            % Number of Decision Variables��Ʊ���ά��
n=nVar;

VarSize=[1 nVar];   % Decision Variables Matrix Size 

lu=[lb;ub];

%% DE Parameters
if n<10
    MaxIt=100;      % Maximum Number of Iterations����������
    nPop=50;        % Population Size��Ⱥ��С
else
    MaxIt=80;      % Maximum Number of Iterations����������
    nPop=100;        % Population Size��Ⱥ��С
end


beta_min=0.2;   % Lower Bound of Scaling Factor  ����F����Сֵ
beta_max=0.8;   % Upper Bound of Scaling Factor  ����F�����ֵ
pCR_min=0.1;    % Lower Bound of Crossover Probability   ����CR����Сֵ 
pCR_max=0.5;    % Lower Bound of Crossover Probability   ����CR����Сֵ



%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

BestSol.Cost=inf;

pop=repmat(empty_individual,nPop,1);

for i=1:nPop

    pop(i).Position=rand()*(lu(2,:)-lu(1,:))+lu(1,:);  %��ʼ����Ⱥ�����λ��������������õĵ�������Ʊ��������½緶Χ�ھ��ȷֲ���˼��
    
    [pop(i).Cost,~]=objfunc(pop(i).Position, par);  %����ÿ����ʼ������Ⱥ�������Ӧ��ֵ
    
    if pop(i).Cost<BestSol.Cost
        BestSol.Cost=pop(i).Cost; %�õ���ʼ����Ⱥ�е���ѵ�λ��������
        BestSol.ind=pop(i).Position;
    end
    
end

BestCost=zeros(MaxIt,1);

%% DE Main Loop

%����DE������
detype=0;%detype=2ʱ�����⣬��Ҫ�Ľ�


FES=nPop;
for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];%�����ǶԵ�i��������б��죬��ˣ�Ҫ��ȥ����i�����壬Ȼ�������
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation  ��ؼ��ı������
        %beta=unifrnd(beta_min,beta_max);
        %����beta�Ƕ���ÿһά���ϵ�ֵ����һ������������ģ�������һ�ε��������и��嶼��һ����ֵ������һ�����������и������Ʊ����Ķ�Ӧά���ϵ�beta��һ����ֵ����ͬά�Ȳ�һ�����Ҳ�ͬ��������Ҳ��һ��
        %��ˣ�����beta��һ��1*Varsize�ľ���
        %�������ò���
        switch detype 
            case 0
                %������򵥵�de,  F��CR������
                beta=0.8;
                pCR=0.8;        % Crossover Probability  �������CR��ֵ
            case 1
                %F��̬�仯����ÿ��ά���Ͼ���ͬ������CR�ǹ̶���
                beta=unifrnd(beta_min,beta_max,VarSize);  %�������������ÿ�ζ�ȡ0.2-0.8��Χ�ڵ�һ�����������ȻҲ���Խ�������Ϊһ����ֵ��beta����F
                pCR=0.2;        % Crossover Probability  �������CR��ֵ
            case 2
                %F��̬�仯����ÿ��ά���Ͼ���ͬ��CR��̬�仯����ÿ��ά���Ͼ���ͬ
                beta=unifrnd(beta_min,beta_max,VarSize);  %�������������ÿ�ζ�ȡ0.2-0.8��Χ�ڵ�һ�����������ȻҲ���Խ�������Ϊһ����ֵ��beta����F
                pCR=unifrnd(pCR_min,pCR_max,VarSize);        % Crossover Probability  �������CR��ֵ
        end
        
        
%         switch mutation
%             case 1
%                 y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);%������õľ��Ǽ򵥵�DE/rand/1�ı��췽����Ҳ�����������ı��췽����Ҫע�Ⲣ������򵥵���ͨ��DE,����������beta��ÿһ��ά�����ǲ�һ����ֵ��
%             case 2
%                 y=BestSol.Position+
%         end
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);%������õľ��Ǽ򵥵�DE/rand/1�ı��췽����Ҳ�����������ı��췽����Ҫע�Ⲣ������򵥵���ͨ��DE,����������beta��ÿһ��ά�����ǲ�һ����ֵ��
        
        VarMin=repmat(lb,size(y,1),1);
        VarMax=repmat(ub,size(y,1),1);
        
        y = max(y, VarMin);
		y = min(y, VarMax);   %ע�����y�Ƕ�Ӧ�ı��������������Ƕ�Ӧ��ĳ��λ�õ�Ŀ�꺯����ֵ���������￼�ǵ������������������ȥ�߽�ȡ�߽磡
        
        
		
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
    
    % Update Best Cost���������Ӧ��ֵ
    BestCost(it)=BestSol.Cost;
    bestvalue=BestSol.Cost;
    bestind=BestSol.ind;
    
    % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end

%% Show Results


end
