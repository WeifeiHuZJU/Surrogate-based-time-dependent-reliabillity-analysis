function  [obj, fit] = objfunc(theta, par)
% t1=clock;
% Initialize
obj = inf; 
fit = struct('sigma2',NaN, 'beta',NaN, 'gamma',NaN, ...
    'C',NaN, 'Ft',NaN, 'G',NaN);
m = size(par.F,1);
% Set up  R
% t10=clock;
r = feval(par.corr, theta, par.D); % 这句代码最耗时，总共评估一次0.06秒的话，这句代码占了0.037秒. 主要原因是 D 矩阵太大。
% t9=clock;
% t_6=etime(t9,t10);
% t_7=etime(t10,t1);
idx = find(r > 0);   o = (1 : m)';   
mu = (10+m)*eps;
% t6=clock;
R = sparse([par.ij(idx,1); o], [par.ij(idx,2); o], ...
  [r(idx); ones(m,1)+mu]);          % 其次是这句代码。总共评估一次0.06秒。这句代码占了0.018秒
% Cholesky factorization with check for pos. def.
% t7=clock;
[C rd] = chol(R);
if  rd,  return, end % not positive definite
% t8=clock;
% Get least squares solution
C = C';   Ft = C \ par.F;
[Q G] = qr(Ft,0);
% t5=clock;
if  rcond(G) < 1e-10
  % Check   F  
  if  cond(par.F) > 1e15 
    T = sprintf('F is too ill conditioned\nPoor combination of regression model and design sites');
    error(T)
  else  % Matrix  Ft  is too ill conditioned
    return 
  end 
end
% t3=clock;
Yt = C \ par.y;   beta = G \ (Q'*Yt);
rho = Yt - Ft*beta;  sigma2 = sum(rho.^2)/m;
detR = prod( full(diag(C)) .^ (2/m) );
obj = sum(sigma2) * detR;
% t2=clock;
% t_Total_time=etime(t2,t1);
% t_latter_time=etime(t2,t3);
% t_former_time=etime(t3,t1);
% t_1=etime(t6,t1);
% t_2=etime(t7,t6);
% t_3=etime(t8,t7);
% t_4=etime(t5,t8);
% t_5=etime(t3,t5);
% save time_record
if  nargout > 1
  fit = struct('sigma2',sigma2, 'beta',beta, 'gamma',rho' / C, ...
    'C',C, 'Ft',Ft, 'G',G');
end