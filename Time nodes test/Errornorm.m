%error norm  误差上确界
% clear all
EN2 = [];
for n = 32:500
%% 先离散化
mu     = 0;        % mean value of the random field   Y（t）的均值
sig    = 1;        % std of the random field  Y（t）的方差
opc    = 1;        % choose example  选择第几个example
switch opc
    case 1 % math example
        corr_leng = 1;    % correlation length x
        cov_func  = @(x) exp(-(x(2)-x(1))^2);
        dom_bound = {[0 1]};
        partition = n;   %  partition in x
end

% solving using EOLE method
[eigval,eigvec,xnod,Mterms] = EOLE_method(cov_func,dom_bound,partition);

% representation of the process: sample function
N         = length(xnod);
NN        = partition;          %
xx        = linspace(dom_bound{1}(1),dom_bound{1}(2),NN);
a = dom_bound{1}(1);
b = dom_bound{1}(2);
parameter.a=a;
parameter.b=b;
parameter.cov_func=cov_func;
cov_mat = zeros(N,NN);
for j = 1:N         %Cy（t）
   for l = 1:NN
      cov_mat(j,l) = cov_func([xnod(j),xx(l)]);
   end
end

%% 计算error norm
mu_u = 0;     %均值
sigma_u = 1;  %方差
mu_u1 = 3.5; 
sigma_u1 = 0.25;

dnum_mcs=1e6;
x1=normrnd(mu_u1,sigma_u1,[dnum_mcs,1]);
x2=normrnd(mu_u1,sigma_u1,[dnum_mcs,1]);
z1=normrnd(mu_u,sigma_u,[dnum_mcs,1]);
z2=normrnd(mu_u,sigma_u,[dnum_mcs,1]);
% z3=normrnd(mu_u,sigma_u,[dnum_mcs,1]);
% t = linspace(a,b,dnum_mcs);
% t = t';
t = unifrnd(a,b,dnum_mcs,1);
z = [z1,z2,t];
parameter.num_mcs1= size(z,1);
Nsim = dnum_mcs;

cov_mat = zeros(Nsim,NN);
for j = 1:Nsim         %Cy（t）
   for l = 1:NN
      cov_mat(j,l) = cov_func([z(j,Mterms+1),xx(l)]);
   end
end

for k = 1:Nsim
    Var=0;
    for i = 1:Mterms
        xi(i,1) = z(k,i);
        Var = Var - (1/eigval(i))*(eigvec(:,i)'*cov_mat(k,:)').^2;
    end
    V_hat(k,1) = sig + Var; 
end
Esuper = mean(V_hat);
EN2 = [EN2,Esuper]
end