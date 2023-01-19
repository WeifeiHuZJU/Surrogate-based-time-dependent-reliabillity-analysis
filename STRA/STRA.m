clear all;clc;
for iii=1:1
    %% Define the problem and initialization
    problem=1;
    N_initial=12;
    total_Nmcs=5e5;
    [X_mean,X_std,x_mcs,z_mcs,Stochasticpro,Xsample]=Case(problem,total_Nmcs,N_initial); % select the type of case study
    Ysample=LSF(problem,Xsample); % true observations
                            
    %% Parameter settings of Kriging moodels
    n_mu=size(Xsample,2);
    theta=repmat(10,1,n_mu);
    lob = repmat(1e-1,1,n_mu);
    upb = repmat(30,1,n_mu); % it will be better to use larger upper bounds, such as 1e3/1e4
%     Using_EA = 0; % Do not use evolutionary algorithm to optimize the hyperparameters of Kriging
    Using_EA = 1; % Use evolutionary algorithm to optimize the hyperparameters of Kriging
    dmodel=dacefit(Xsample, Ysample, @regpoly1, @corrgauss, theta, lob, upb,Using_EA);
    
    %% identify some training points
    p=0;error_max_set=[];error_min_set=[];error_max=1;pf_set=[];error_target=0.05;N_set=[];Prob_wrong_max=[];    
    while p==0 || error_max>error_target
        p=p+1
        dmodel=dacefit(Xsample, Ysample, @regpoly1, @corrgauss, theta,lob,upb,Using_EA); % update Kriging model
        theta_X=dmodel.theta; % correlation function parameters
        % Calculate U function for each random point
        U1=[];U2=[];g_x_t1=[];g_x_t2=[];Ft=[];MCS_t=[];MCS_t_normal=[];X_normal=[];
        for t=1:Stochasticpro.Nnode
            tnode=Stochasticpro.tnode;
            T=Stochasticpro.T;
            Ft=EOLE_out(tnode(t),z_mcs,Stochasticpro);
            if size(Xsample,2)-length(X_mean)==1
                MCS_t=[x_mcs,Ft];
            elseif size(Xsample,2)-length(X_mean)==2
                MCS_t=[x_mcs,Ft,repmat(tnode(t),total_Nmcs,1)];
            end
            for i=1:(total_Nmcs/(2*1e4))
                [g_x_t1((1+(i-1)*1e4):i*1e4,t),mse_x_t1((1+(i-1)*1e4):i*1e4,t)]=predictor(MCS_t((1+(i-1)*1e4):i*1e4,:),dmodel);
                [g_x_t2((1+(i-1)*1e4):i*1e4,t),mse_x_t2((1+(i-1)*1e4):i*1e4,t)]=predictor(MCS_t((total_Nmcs/2+1+(i-1)*1e4):(total_Nmcs/2+i*1e4),:),dmodel);
            end
            U1(:,t)=abs(g_x_t1(:,t))./sqrt(mse_x_t1(:,t));U2(:,t)=abs(g_x_t2(:,t))./sqrt(mse_x_t2(:,t));
            
            %% 相关性标准核验 compute the correlation
            if size(Xsample,2)-length(X_mean)==1
                Wmax=[max(x_mcs),max(Ft)];
                Wmin=[min(x_mcs),min(Ft)];
            elseif size(Xsample,2)-length(X_mean)==2
                Wmax=[max(x_mcs),max(Ft),T(2)];
                Wmin=[min(x_mcs),min(Ft),T(1)];
            end
            MCS_t_normal=2*(MCS_t-repmat(Wmin,total_Nmcs,1))./repmat(Wmax-Wmin,total_Nmcs,1)-1;
            X_normal=2*(Xsample-repmat(Wmin,size(Xsample,1),1))./repmat(Wmax-Wmin,size(Xsample,1),1)-1;
            psi=[];psi_max=[];psi_t1=[];psi_t2=[];
            for j=1:size(Xsample,1)
                psi(:,j)=exp(-sum(repmat(theta_X,total_Nmcs,1).*abs(repmat(X_normal(j,:),total_Nmcs,1)-MCS_t_normal).^2,2));
             end
            psi_max=max(psi');%psi每行的最大值
            psi_t1=find(psi_max(1:total_Nmcs/2)>=0.99);psi_t2=find(psi_max(total_Nmcs/2+1:total_Nmcs)>=0.99);
            U1(psi_t1,t)=100;U2(psi_t2,t)=100; %剔除不满足相关性要求的点
        end
        
        %% 计算失效概率
        nf=0;
        min_g_x_t1=min(g_x_t1');min_g_x_t2=min(g_x_t2');
        min_g_x_t1=min_g_x_t1';min_g_x_t2=min_g_x_t2';
        min_g_x_t=[min_g_x_t1;min_g_x_t2];
        nf=length(find(min_g_x_t<=0));
        pf=nf/total_Nmcs;
        pf_set=[pf_set;pf]
        
        %% 计算失效概率误差
        U_t1=U1;U_t2=U2;
%         U_t1=abs(g_x_t1)./sqrt(mse_x_t1);U_t2=abs(g_x_t2)./sqrt(mse_x_t2);%U_t1和Ut2是没有加相关性约束的U
        Prob_wrong_t1=normcdf(-U_t1);Prob_wrong_t2=normcdf(-U_t2);
        ig=ones(total_Nmcs,1);
        N_s_predict=sum(ig(min_g_x_t>0));
        Prob_wrong_t_s1=Prob_wrong_t1(min_g_x_t1>0,:);Prob_wrong_t_s2=Prob_wrong_t2(min_g_x_t2>0,:);%把候选点分为安全点和失效点
        Prob_right_t_s1=1-Prob_wrong_t_s1;Prob_right_t_s2=1-Prob_wrong_t_s2;
        Prob_wrong_s1=1-prod(Prob_right_t_s1,2);Prob_wrong_s2=1-prod(Prob_right_t_s2,2);
        Prob_right_s1=prod(Prob_right_t_s1,2);Prob_right_s2=prod(Prob_right_t_s2,2);
        Prob_wrong_s=[Prob_wrong_s1;Prob_wrong_s2];
        Prob_right_s=[Prob_right_s1;Prob_right_s2];
        mean_wrong_N_s=sum(Prob_wrong_s);%Kriging预测的安全域内错误分类的候选点个数的均值
        std_wrong_N_s=sqrt(sum(Prob_wrong_s.*Prob_right_s));%Kriging预测的安全域内错误分类的候选点个数的标准差
        %计算错分的失效点的均值和方差
        N_f_predict=sum(ig(min_g_x_t<=0));
        Prob_wrong_t_f1=Prob_wrong_t1(min_g_x_t1<=0,:);Prob_wrong_t_f2=Prob_wrong_t2(min_g_x_t2<=0,:);
        Prob_wrong_t_f=[Prob_wrong_t_f1;Prob_wrong_t_f2];
        Prob_right_t_f=1-Prob_wrong_t_f;
        g_x_t_f=[g_x_t1(min_g_x_t1<=0,:);g_x_t2(min_g_x_t2<=0,:)];%失效点每个时刻的值
        Leicheng1=[];Leicheng2=[];
        for ii=1:N_f_predict
            Leicheng1(ii,1)=prod(Prob_wrong_t_f(ii,g_x_t_f(ii,:)<=0));
            Leicheng2(ii,1)=prod(Prob_right_t_f(ii,g_x_t_f(ii,:)>0));
        end       
        Prob_wrong_f=Leicheng1.*Leicheng2;
        Prob_right_f=1-Prob_wrong_f;
        mean_wrong_N_f=sum(Prob_wrong_f);%Kriging预测的失效域内错误分类的候选点个数的均值
        std_wrong_N_f=sqrt(sum(Prob_wrong_f.*Prob_right_f));%Kriging预测的失效域内错误分类的候选点个数的标准差
        
        N_s_up=norminv(0.975,mean_wrong_N_s,std_wrong_N_s); %Kriging预测的安全域内错误分类的候选点个数的上界
        N_s_low=norminv(0.025,mean_wrong_N_s,std_wrong_N_s);
        N_f_up=poissinv(0.975,mean_wrong_N_f); %Kriging预测的失效域内错误分类的候选点个数的上界
        N_f_low=poissinv(0.025,mean_wrong_N_f);
        N_s_f=[N_s_low,N_s_up,N_f_low,N_f_up];
        N_set=[N_set;N_s_f];
        % 两种计算误差的方式，只有微小差别
%         error_max=max(abs(N_f_predict./(N_f_predict+N_s_low-N_f_up)-1),abs(N_f_predict./(N_f_predict+N_s_up-N_f_low)-1));
        error_max=max(abs(N_f_predict./(N_f_predict-N_f_up)-1),abs(N_f_predict./(N_f_predict+N_s_up)-1));
        error_max_set=[error_max_set;error_max]
        if error_max<error_target
            break;
        end
        
        %% 确定新的样本点
        t = linspace(T(1),T(2),total_Nmcs);
        Ft_loo=EOLE_out(t,z_mcs,Stochasticpro);
        sam_domain_x = [x_mcs,Ft_loo,t'];
        [cellposition,weight]=LOO_5(Xsample, Ysample,dmodel,sam_domain_x,Stochasticpro);
        sample=[];
        mWEF=[];
         for i = 1:size(cellposition,1)
            %划分voronoi cell，并且筛选出在最敏感分区中的候选样本点

            [cell]=voronoi(sam_domain_x,Xsample,cellposition(i));
             if (size(cell,1)==0||(size(cell,1)==1))
                continue
            end
            [y_kr, mse_kr] = predictor(cell, dmodel);
            EF = EFF(y_kr, abs(mse_kr));
            WEF= EF* weight(i);
            mWEF(i) = max(WEF);
            sam_center = find(WEF==max(WEF));
            sample = [sample;cell(sam_center,:)];
         end
         sam_location = find(mWEF==max(mWEF));
         X_new = sample(sam_location,:);
  
        Y_new=LSF(problem,X_new);
        Xsample=[Xsample;X_new];
        Ysample=[Ysample;Y_new];
  
    end
    % 计算最终结果
    Nnode=Stochasticpro.Nnode;
    nf_t=zeros(Nnode,1);
    ig1=ones(total_Nmcs/2,1);
    ig2=ones(total_Nmcs/2,1);
    for i=1:Nnode
        nf1(i,1)=sum(ig1(min(g_x_t1(:,1:i),[],2)<0));
        nf2(i,1)=sum(ig2(min(g_x_t2(:,1:i),[],2)<0));
    end
    nf_t=nf1+nf2;
    pf_t=nf_t/total_Nmcs;%失效概率函数
    min_g_x_t1=min(g_x_t1');min_g_x_t2=min(g_x_t2');
    min_g_x_t_all(:,iii)=[min_g_x_t1';min_g_x_t2'];
    Number_X(iii)=size(Xsample,1);
    pfff(:,iii)=pf_t;
    
    for i=1:Stochasticpro.Nnode
        tnode=Stochasticpro.tnode;
        Ft=EOLE_out(tnode(i),z_mcs,Stochasticpro);
        if size(Xsample,2)-length(X_mean)==1
            MCS_t=[x_mcs,Ft];
        elseif size(Xsample,2)-length(X_mean)==2
            MCS_t=[x_mcs,Ft,repmat(tnode(i),total_Nmcs,1)];
        end
        g1_true(:,i)=LSF(problem,MCS_t(1:total_Nmcs/2,:));
        g2_true(:,i)=LSF(problem,MCS_t(total_Nmcs/2+1:end,:));
    end
    nf=0;
    min_g_x_t1_true=min(g1_true');
    min_g_x_t2_true=min(g2_true');
    nf_true=length(find(min_g_x_t1_true<0))+length(find(min_g_x_t2_true<0));
    pf_true=nf_true/total_Nmcs;
    error_true=abs(pf_set-pf_true)./pf_true;
    Error{iii}=[error_max_set,error_true];
    
end