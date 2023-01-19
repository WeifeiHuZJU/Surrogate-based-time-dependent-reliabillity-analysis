function [cellposition,weight]=LOO_5(sam_x,sam_y,model_cons1,sam_domain_x,Stochasticpro)
%leave one out technique

m = size(sam_x,1);%number of doe
sam_domain_x_loo = sam_domain_x;
NN = Stochasticpro.Nnode;
error=zeros(m,1);

% [sam_domain_y]=cons2(sam_domain_x,parameter);%采样点对应的响应值
[y_kr, ~] = predictor(sam_domain_x, model_cons1);
I = zeros(length(sam_domain_x_loo),NN);
I(y_kr < 0) = 1;
I_p = sum(sum(I,2)>0);
p_f = I_p/(length(sam_domain_x_loo));

% rmse_all = sqrt(mean((sam_domain_y - y_kr).^2));

% sam_domain_y_loo = sam_domain_y;

%LOO
for i=1:m
    sam_x_loo=sam_x;
    sam_y_loo=sam_y;
    sam_x_loo(i,:)=[];%去除第i个样本
    sam_y_loo(i,:)=[];
    parameter.num_mcs1= size(sam_x_loo,1);
    parameter.num_doe_cons1 = size(sam_x_loo,1);

     n_mu=size(sam_x,2);
    theta=repmat(10,1,n_mu);
    lob = repmat(1e-1,1,n_mu);
    upb = repmat(30,1,n_mu); % it will be better to use larger upper bounds, such as 1e3/1e4
%     Using_EA = 0; % Do not use evolutionary algorithm to optimize the hyperparameters of Kriging
    Using_EA = 1; % Use evolutionary algorithm to optimize the hyperparameters of Kriging
    model_LOO=dacefit(sam_x_loo,sam_y_loo, @regpoly1, @corrgauss, theta, lob, upb,Using_EA);
    
    [y_kr_LOO, ~] = predictor(sam_domain_x_loo, model_LOO);
    I = zeros(length(sam_domain_x_loo),NN);
    I(y_kr_LOO < 0) = 1;
    I_p_LOO = sum(sum(I,2)>0);
    p_f_LOO = I_p_LOO/(length(sam_domain_x_loo));
    
    error(i)=abs(p_f_LOO-p_f);
end

% [AS,pos]=sort(error(:),'descend');
% cellposition = pos;

mean_error = mean(error);
[errormax,~]=max(error);

[AS,pos]=sort(error(:),'descend');
for j = 1:size(error,1)
    if AS(j,1)> mean_error
        cellposition(j,1) = pos(j,1);
        weight(j,1)=AS(j)/errormax;
    end
end



end
