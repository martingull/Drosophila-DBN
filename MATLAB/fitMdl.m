function mdl = fitMdl(X, y, para) 
% Computes coefficients for network. 

if nargin < 3
    para = get_parameters() ;
end

if strcmpi(para.robustOpts, 'gaussian') % Gaussian ML
    mdl = fitlm(X, y) ; 
    %mdl = linregFit(X, y) ; % same source as linregFitBayes
    
elseif strcmpi(para.robustOpts, 'EB') % Empirical Bayes
    [mdl, logev] = linregFitBayes(X, y, 'prior', 'eb') ;
    mdl.LogLikelihood = logev ;
    
elseif strcmpi(para.robustOpts, 'VB') % Variational Bayes
    [mdl, logev] = linregFitBayes(X, y, 'prior', 'vb') ;
    mdl.LogLikelihood = logev ;
    
elseif strcmpi(para.robustOpts, 'constML') % Constrained Maximum Likelihood
    mdl = fitConstML(X, y) ;
    %mdl.LogLikelihood = logev ;

else
    disp('Incorrect robustOpts passed from get_parameters.')    
end


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

