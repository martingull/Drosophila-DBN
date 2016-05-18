function [predData, score] = mk_simulation(data, inter, para)
% Simulate data given data and causal structure.
% data(gene,instance,time) 

% If para not given then global parameters will be used.
% OBS: data_stack used global parameters. 
if nargin < 3
    para = get_parameters() ;
end

% Extract Parameters
endo = para.endo ;
[G, N, T] = size(data);
nEndo = length(endo) ;


[dataFix, yRHS, yLHS] = mk_dataStack(data, para) ;


% Training using intertermporal graph
mdl = cell(1, nEndo) ;
parents = cell(1, max(endo)) ;
RMSE = zeros(1, nEndo) ;
score = zeros(1, nEndo) ;


for i = endo
    % Time Parameters
    
    parents{i} = find(inter(:,i))' ;
    X = [dataFix{i} yRHS(:, parents{i})] ;
    
    % fitlm(X, y) ~ y = bX + eps
    mdl{i} = fitMdl(X, yLHS(:, i), para) ; 

    % Compute simulation noise.
    RMSE(i) = mk_RMSE(mdl{i}, X, para) ;
    
    score(i) = mk_score(mdl{i}, para) ;
end

score = sum(score) ;


% Predict, all data, overwrite endo in the simulation.
predData = nan(G,N,T,para.nSim) ;

% parfor i = 1:para.nSim
for i = 1:para.nSim
simData = data ; % (gene x space x time)
    for t = 2:T
        for g = endo 

            %X = [simData(g,2:end-1,t-1)' mk_autoregCube(simData(g,2:end-1,t-1)',para) mk_euler(simData(g,:,t-1)',para) simData(parents{g},2:end-1,t-1)'] ;
            X = [simData(g,:,t-1)' mk_euler(simData(g,:,t-1)',para) simData(parents{g},:,t-1)'] ;
            %X = [simData(g,2:end-1,t-1)' mk_euler(simData(g,:,t-1)',para) simData(parents{g},2:end-1,t-1)'] ;
            %X = [simData(g,:,t-1)' mk_autoregCube(simData(g,:,t-1)',para) mk_euler(simData(g,:,t-1)',para) simData(parents{g},:,t-1)'] ;

            % predict: y ~ Xb + eps
            simData(g,:,t) = predMdl(mdl{g}, X, RMSE(g), para) ;

        end
    end
    
predData(:,:,:,i) = simData ;
end
    

end


function yPred = predMdl(mdl, X, noise, para)
% Predicts using pred function.

if nargin < 4
    para = get_parameters() ;
end

if strcmpi(para.robustOpts, 'gaussian')    
	yPred = random(mdl, X); % add noise
    
elseif strcmpi(para.robustOpts, 'EB') || strcmpi(para.robustOpts, 'VB')
    [yPred, noise] = linregPredictBayes(mdl, X) ;
    yPred = normrnd(yPred, noise) ; % add noise
    
elseif    strcmpi(para.robustOpts, 'constML')
    yPred = normrnd([ones(size(X,1),1) X] * mdl.theta, noise) ;
    
else
    disp('predMdl() passes nonsupported robustOpts.')
    
end


end

function RMSE = mk_RMSE(mdl, X, para) 
% Computes Root mean squared error for mdl objects.
if nargin < 3
    para = get_parameters() ;
end


if strcmpi(para.robustOpts, 'gaussian')
    RMSE = mdl.RMSE ;
   
elseif strcmpi(para.robustOpts, 'EB') 
    RMSE = sqrt(mdl.sigma2) ;
    
elseif strcmpi(para.robustOpts, 'VB')
    % Use proxy function to compute.
    [~, sig2] = linregPredictBayes(mdl, X) ;
    RMSE = sqrt(mean(sig2)) ;
    
elseif    strcmpi(para.robustOpts, 'constML')
    RMSE = mdl.sigma ;
    
    
end
    
    
end

function score = mk_score(mdl, para)

if strcmpi(para.robustOpts, 'gaussian')
    score = mdl.criterion.BIC ;
   
elseif strcmpi(para.robustOpts, 'EB') 
    score = mdl.LogLikelihood ; % log-evidence
    
elseif strcmpi(para.robustOpts, 'VB')
    score = mdl.LogLikelihood ; % log-evidence
    
elseif    strcmpi(para.robustOpts, 'constML')
    score = mdl.BIC ;
        
end

end

