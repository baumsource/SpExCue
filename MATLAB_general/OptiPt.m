function [p,chistat,u,lL_eba,lL_sat,fit,cova] = OptiPt(M,A,s)
% OptiPt parameter estimation for BTL/Pretree/EBA models
%   p = OptiPt(M,A) estimates the parameters of a model specified
%   in A for the paired-comparison matrix M. M is a matrix with
%   absolute frequencies. A is a cell array.
%
%   [p,chistat,u] = OptiPt(M,A) estimates parameters and reports
%   the chi2 statistic as a measure of goodness of fit. The vector
%   of scale values is stored in u.
%
%   [p,chistat,u,lL_eba,lL_sat,fit,cova] = OptiPt(M,A,s) estimates
%   parameters, checks the goodness of fit, computes the scale values,
%   reports the log-likelihoods of the model specified in A and of the
%   saturated model, returns the fitted values and the covariance
%   matrix of the parameter estimates. If defined, s is the starting
%   vector for the estimation procedure. Otherwise each starting value
%   is set to 1/length(p).
%   The minimization algorithm used is FMINSEARCH.
%
%   Examples
%     Given the matrix M = 
%                            0    36    35    44    25
%                           19     0    31    37    20
%                           20    24     0    46    24
%                           11    18     9     0    13
%                           30    35    31    42     0
%
%     A BTL model is specified by A = {[1];[2];[3];[4];[5]}
%     Parameter estimates and the chi2 statistic are obtained by
%       [p,chistat] = OptiPt(M,A)
%
%     A Pretree model is specified by A = {[1 6];[2 6];[3 7];[4 7];[5]} 
%     A starting vector is defined by s = [2 2 3 4 4 .5 .5]
%     Parameter estimates, the chi2 statistic, the scale values, the
%     log-likelihoods of the Pretree model and of the saturated model,
%     the fitted values, and the covariance matrix are obtained by
%       [p,chistat,u,lL_eba,lL_sat,fit,cova] = OptiPt(M,A,s)
%
% Authors: Florian Wickelmaier (wickelmaier@web.de) and Sylvain Choisel
% Last mod: 03/JUL/2003
% For detailed information see Wickelmaier, F. & Schmid, C. (2004). A Matlab
% function to estimate choice model parameters from paired-comparison data.
% Behavior Research Methods, Instruments, and Computers, 36(1), 29-40.

I = length(M);  % number of stimuli
mmm = 0;
for i = 1:I
  mmm = [mmm max(A{i})];
end
J = max(mmm);  % number of pt parameters
if(nargin == 2)
  p = ones(1,J)*(1/J);  % starting values
elseif(nargin == 3)
  p = s;
end

for i = 1:I
  for j = 1:I
    diff{i,j} = setdiff(A{i},A{j});  % set difference
  end
end

p = fminsearch(@ebalik,p,optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000),M,diff,I);  % optimized parameters
lL_eba = -ebalik(p,M,diff,I);  % likelihood of the specified model

lL_sat = 0;  % likelihood of the saturated model
for i = 1:I-1
  for j = i+1:I
    lL_sat = lL_sat + M(i,j)*log(M(i,j)/(M(i,j)+M(j,i)))...
                    + M(j,i)*log(M(j,i)/(M(i,j)+M(j,i)));
  end
end

fit = zeros(I);  % fitted PCM
for i = 1:I-1
  for j = i+1:I
    fit(i,j) = (M(i,j)+M(j,i))/(1+sum(p(diff{j,i}))/sum(p(diff{i,j})));
    fit(j,i) = (M(i,j)+M(j,i))/(1+sum(p(diff{i,j}))/sum(p(diff{j,i})));
  end
end

chi = 2*(lL_sat-lL_eba);
df =  I*(I-1)/2 - (J-1);
chistat = [chi df];  % 1-chi2cdf(chi,df)];  % goodness-of-fit statistic

u = sum(p(A{1}));  % scale values
for i = 2:I
  u = [u sum(p(A{i}))];
end

H = hessian('ebalik',p',M,diff,I);
C = inv([H ones(J,1); ones(1,J) 0]);
cova = C(1:J,1:J);
end

function lL_eba = ebalik(p,M,diff,I)  % computes the likelihood

if min(p)<=0  % bound search space
  lL_eba = inf;
  return
end

thesum = 0;
for i = 1:I-1
  for j = i+1:I
    thesum = thesum + M(i,j)*log(1+sum(p(diff{j,i}))/sum(p(diff{i,j})))...
                    + M(j,i)*log(1+sum(p(diff{i,j}))/sum(p(diff{j,i})));
  end
end
lL_eba = thesum;
end

function H = hessian(f,x,varargin)  % computes numerical Hessian
% based on a solution posted on Matlab Central by Paul L. Fackler

k = size(x,1);
fx = feval(f,x,varargin{:});
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;
ee = sparse(1:k,1:k,h,k,k);

g = zeros(k,1);
for i = 1:k
  g(i) = feval(f,x+ee(:,i),varargin{:});
end

H = h*h';
for i = 1:k
  for j = i:k
    H(i,j) = (feval(f,x+ee(:,i)+ee(:,j),varargin{:})-g(i)-g(j)+fx)...
                 / H(i,j);
    H(j,i) = H(i,j);
  end
end
end