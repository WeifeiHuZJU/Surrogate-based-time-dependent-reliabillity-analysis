function  [r, dr] = corrgauss(theta, d)
%CORRGAUSS  Gaussian correlation function,
%
%           n
%   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = corrgauss(theta, d)
%          [r, dr] = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update June 2, 2002

[m n] = size(d);  % number of differences and dimension of data  这里若 样本点数为 mm 。 则  m=mm*(mm-1)/2 。所以样本点数增加，m增加巨大。
if  length(theta) == 1
  theta = repmat(theta,1,n);
elseif  length(theta) ~= n
  error(sprintf('Length of theta must be 1 or %d',n))
end

td = d.^2 .* repmat(-theta(:).',m,1);   % 这句代码就是 Kriging 建模最耗时的那句代码！！！！！

r = exp(sum(td, 2));

if  nargout > 1
  dr = repmat(-2*theta(:).',m,1) .* d .* repmat(r,1,n);
end



