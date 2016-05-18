function mdl = fitMdl(X, y, para) 
% Computes coefficients for network. 

if nargin < 3
    para = get_parameters() ;
end

if strcmpi(para.robustOpts, 'gaussian') % Gaussian ML
    mdl = fitlm(X, y) ; 
    %mdl = linregFit(X, y) ; % same source as linregFitBayes
    
elseif strcmpi(para.robustOpts, 'EB') % Empirical Bayes
    mdl = linregFitEbChen(X, y) ;
    mdl.LogLikelihood = logev ;
    
elseif strcmpi(para.robustOpts, 'VB') % Variational Bayes
    % [mdl, logev] = linregFitBayes(X, y, 'prior', 'vb') ;
    % mdl.LogLikelihood = logev ;
    mdl = vb_linear_fit(X, y) ; % Add mdl struct for method 
     
elseif strcmpi(para.robustOpts, 'constML') % Constrained Maximum Likelihood
    mdl = fitConstML(X, y) ;
    %mdl.LogLikelihood = logev ;

else
    disp('Incorrect robustOpts passed from get_parameters.')    
end


end


function [model, L, Lhist] = linregFitEbChen(X, Y)
% Evidence procedure (Empirical Bayes) for linear regression
% PMTKauthor Tao Chen
% PMTKurl http://www3.ntu.edu.sg/home/chentao/software.htm
% It estimates alpha (scalar) and beta
% L is log marginal likelihood
% gamma is effective number of paramers

% This file is from pmtk3.googlecode.com


% we currently ignore whether we prepended 1s  to X or not

[N, M] = size(X);
X = [ones(N, 1) X] ;
M = M + 1 ; 


% pre-computation & initial values

XX = X'*X;
XX2 = X*X';
Xy = X' * Y;

% The method can get stuck in local minima, so we should
% do multiple restarts
%alpha = exp(randn()*3-3); %alpha=1;
%beta = exp(randn()*3-3); %beta=1;
alpha = 0.01; % initiailly don't trust prior
beta = 1; % initially trust the data
mn = zeros(M,1); Sn = zeros(M,M);

L_old = -inf;
Lhist = [];
for i = 1:100;
  
  % calcualte covariance matrix S
  if ( N > M )
    T = alpha*eye(M) + XX*beta;
    cholT = chol(T);
    Ui = inv(cholT);
    Sn = Ui * Ui';
    logdetS = - 2 * sum ( log(diag(cholT)) );
  else
    T = eye(N)/beta + XX2/alpha;
    cholT = chol(T);
    Ui = inv(cholT);
    Sn = eye(M)/alpha - X' * Ui * Ui' * X / alpha / alpha;
    logdetS = sum(log(diag(cholT)))*2 + M*log(alpha) + N*log(beta);
    logdetS = - logdetS;
  end
  
  mn = beta * Sn * Xy;
 
  
  t1 = sum ( (Y - X * mn).^2 );
  t2 = mn' * mn;
  
  gamma = M - alpha * trace(Sn);
  beta = ( N - gamma ) / ( t1 );
  
  L = M*log(alpha) - N*log(2*pi) + N*log(beta) - beta*t1 - alpha*t2 + logdetS;
  L = L/2;
  Lhist(i) = L; %ok
  fprintf('Iter %d: L=%f, alpha=%f, beta=%f\n', i, L, alpha, beta);
  
  
  if abs(L - L_old) < 1e-2 % use absolute change to avoid small uphill steps
      if (mn(2)-4*mn(3))^2 < 1 && mn(2)^2 < 1
          break;
      else
          alpha = 1.2 * alpha ;
      end%  especially at the initial iterations
  else
    % update alpha only if we DO NOT break
    alpha = ( gamma ) / ( t2 );
  end
  L_old = L;
  
  
  
end
% Needed by predict
model.wN = mn;
model.VN = Sn;
model.beta = beta;

% For diagnostic purposes only
model.alpha = alpha;
model.gamma = gamma;


end






function [w, V, invV, logdetV, an, bn, E_a, L] = vb_linear_fit(X, y, a0, b0, c0, d0)
% [w, V, invV, logdetV, an, bn, E_a, L] = vb_linear_fit(X, y)
%
% estimates w sucht that y = Xw, using Bayesian regularisation.
%
% The underlying generative model assumes
%
% p(y | x, w, tau) = N(y | w'x, tau^-1),
%
% with x and y being the rows of the given X and y. w and tau are assigned
% the conjugate normal inverse-gamma prior
%
% p(w, tau | alpha) = N(w | 0, (tau alpha)^-1 I) Gam(tau | a0, b0),
%
% with the hyper-prior
%
% p(alpha) = p(alpha | c0, d0).
%
%
% The prior parameters a0, b0, c0, and d0 can be set by calling the script
% with the additional parameters vb_linear_fit(X, y, a0, b0, c0, d0). If
% not given, they default to values a0 = 1e-2, b0 = 1e-4, c0 = 1e-2, and
% d0 = 1e-4, such that the prior is uninformative.
%
% The returned posterior parameters (computed by variational Bayesian
% inference) determine a posterior of the form
%
% N(w1 | w, tau^-1 V) Gam(tau | an, bn).
%
% Also, the mean E_a = E(alpha) is returned, together with the inverse of V,
% and its log determinant. L is the variational bound of the model, and is a
% lower bound on the log-model evidence ln p(y | X).
%
% Copyright (c) 2013, 2014, Jan Drugowitsch
% All rights reserved.
% See the file LICENSE for licensing information.



% prior parameters Uninformative
if nargin < 3,  a0 = 1e-6;  end
if nargin < 4,  b0 = 1e-6;  end
if nargin < 5,  c0 = 1e-6;  end
if nargin < 6,  d0 = 1e-6;  end


% pre-process data
[N, D] = size(X);
X = [ones(N, 1) X] ; 
D = D + 1;
X_corr = X' * X;
Xy_corr = X' * y;
an = a0 + N / 2;    gammaln_an = gammaln(an);
cn = c0 + D / 2;    gammaln_cn = gammaln(cn);

% iterate to find hyperparameters
L_last = -realmax;
max_iter = 500;
E_a = c0 / d0;
for iter = 1:max_iter
    % covariance and weight of linear model
    invV = E_a * eye(D) + X_corr;
    V = inv(invV);
    logdetV = - logdet(invV);
    w = V * Xy_corr;
    % parameters of noise model (an remains constant)
    sse = sum((X * w - y) .^ 2);
    bn = b0 + 0.5 * (sse + E_a * (w' * w));
    E_t = an / bn;

    % hyperparameters of covariance prior (cn remains constant)
    dn = d0 + 0.5 * (E_t * (w' * w) + trace(V));
    E_a = cn / dn;

    % variational bound, ignoring constant terms for now
    L = - 0.5 * (E_t * sse + sum(sum(X .* (X * V)))) + 0.5 * logdetV ...
        - b0 * E_t + gammaln_an - an * log(bn) + an ...
        + gammaln_cn - cn * log(dn);

    % variational bound must grow!
    if L_last > L
        fprintf('Last bound %6.6f, current bound %6.6f\n', L_last, L);
        error('Variational bound should not reduce');
    end
    % stop if change in variation bound is < 0.001%
    if abs(L_last - L) < abs(0.00001 * L) 
        break
    end
    L_last = L;    
end
if iter == max_iter
    warning('Bayes:maxIter', ...
        'Bayesian linear regression reached maximum number of iterations.');
end


% augment variational bound with constant terms
L = L - 0.5 * (N * log(2 * pi) - D) - gammaln(a0) + a0 * log(b0) ...
    - gammaln(c0) + c0 * log(d0);
end


function mdl = fitConstML(X, y)
% Compute constrained ML (eigenvalues within unit circle)

[nObs, k] = size(X) ;
X = [ones(nObs, 1) X] ; 
k = k + 1;


% Set Options for fmincon, TolX is the stepsize tollerance
options = optimoptions('fmincon', 'Display', 'off','TolX',1.0e-16);%,'GradObj','on');


% Define objective function -log(likelihood)
f =@(theta)objFunc(theta, X, y) ;
%fb =@(theta)gradient(theta, X, y) ; 


% Initial conditions of search(cheap regression)
x0 = regress(y, X)' ; %

% Constraints %%%%%%%%%%%%%%%%%%%%

% For nhood diffusion 
%A = zeros(1,k) ;
%A(:,2) = 1 ;
%A(1,3) = 2 ;
%A(2,3) = -2 ;
%b = [0.99; 0.99] ;

% For euler diffusion
A = [] ;
b = [] ;
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon =@ eigStab;
[theta, fval] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon,options) ;
% fval ~ smallest negative log-likelihood  
% x ~ parameters
%%%%%%%%%%%%%%

% Transform variables, u_new = exp(u_old) i.e u_new>0
% If you're lost Google: Box Cox Transformations.
%theta(2) = exp(theta(2)) ; % 1 > AR term > -1
theta(3) = exp(theta(3)) ; % diffusions term > 0

mdl.theta = theta' ;
mdl.LogLikelihood = -fval ;
mdl.BIC = -2*mdl.LogLikelihood + k.*(log(nObs)-log(2*pi)) ;
mdl.sigma = std(y - X*theta') ; % standard deviation of residuals.


end

function [c,ceq] = eigStab(theta)
%c = (exp(theta(2)) - 4 * exp(theta(3)))^2 - 0.99999 ;
c = (theta(2) - 4 * exp(theta(3)))^2 - 0.99999 ;
ceq = [];
end

function [f] = objFunc(theta, X, y)
    % Objective function log-likelihood.
    [n, ~] = size(X) ;
    
    % Precomputation. 
    %mu = X*[theta(1) exp(theta(2)) exp(theta(3)) theta(4:end)]'  ;
    mu = X*[theta(1) theta(2) exp(theta(3)) theta(4:end)]'  ;
    z = (y - mu) ;
    
    % Use ML estimator for sigma
    sigma = sqrt((z' * z) / (n - 1)) ;
    % Compute log-likelihood.
    f = -sum(log(normpdf(y, mu, sigma))) ;
    
    % Compute gradient.
    %g = (z' * X) / (sigma * sigma) ;
    % to add gradient (accuracy/speedup) make argument out [f, g]
    % 


end

