function inter = search_dbn(logData, para) 
% This function learns the structure of the intertemporal DBN. 



% Get global parameters
if nargin < 2
    para = get_parameters() ;
end
    
% Predim
interCell = cell(1, para.nRestart);
BIC = zeros(1, para.nRestart);

% CAN ADD PARFOR HERE TO USE MORE CORES. 
parfor i = 1:para.nRestart
     [interCell{i}, BIC(i)] = learn_struct_dbn_stats(logData, para) ;
     disp(['Run number ', num2str(i), ' of ', num2str(para.nRestart), ' complete.'])
end

% Find smallest BIC structure.
[~, bestBIC] = min(BIC) ; % max LogL <-> min BIC
inter = interCell{bestBIC} ;
end

function [inter, BIC] = learn_struct_dbn_stats(data, para) 


% Get global parameters
if nargin < 2
    para = get_parameters() ;
end
endo = para.endo ;


% Initial sizes
[G, N, T] = size(data);
nEndo = length(endo); 
nObs = (N-2) * (T-1) ;

% Initialise Inter
inter = zeros(G);
inter(:, para.endo) = randi([0,1], G, length(para.endo)) ; % random init
inter(logical(eye(G))) = 0 ; % selfloop included in fixed effects.


% Reformulate data 
[dataFix, yRHS, yLHS] = mk_dataStack(data) ;


% Compute spatial likelihood


% compute number of fixed parameters
cnt = nEndo * (size(dataFix{1}, 2) + 2) ; % nEndo(c + ar + causes).


% Initial Scores (4 endo scores), only diff/feedforward
mdl = cell(1, G) ;
logL = zeros(1, G) ;
for i = endo
    mdl{i} = fitMdl(dataFix{i}, yLHS(:, i), para) ; % fitlm(X, y) ~ y = bX + eps
    logL(i) = mdl{i}.LogLikelihood ; 
end
%[~,BIC] = aicbic(sum(logL), cnt, nObs) ;
BIC = mdlScore(sum(logL), cnt, nObs) ;


% GREED ALGO 

for i = 1:para.nIter
    
    % Pick endo node
    node = endo(randi([1, nEndo])); 
    
    % Find Parents, exclude oneself and current parents
    parents = find(inter(:, node))' ;
    avail_parents = setdiff(setdiff(1:G , parents), node) ; 
    
    
    if rand() > 0.5  % Add Parent 

        % Find lowest BIC-score
        jLL = -Inf(1, G) ;
        for j = avail_parents %intersect(1:nEndo, avail_parents)
             % X = Diffusion/Feedforward + current parents + suggested parent
             X = [dataFix{node} yRHS(:, parents) yRHS(:, j) ];
             jmdl = fitMdl(X, yLHS(:, node)) ;
             jLL(j) = jmdl.LogLikelihood ; 
        end
        [iLL, iParent] = max(jLL) ; 

        % Creat new logL vector
        ilogL = logL ;
        ilogL(node) = iLL ;
        nParams = cnt + sum(sum(inter)) + 1 ;
        iBIC = mdlScore(sum(ilogL), nParams, nObs) ;
        
        % See if BIC is lower, if lower approve
        if iBIC < BIC 
            inter(iParent, node) = 1 ; % inter(PARENT, NODE)
            BIC = iBIC ;
            logL = ilogL ;
            %disp(['Gene ', num2str(node), ' got parent ', num2str(iParent)])
        end
        
        
    else % Remove Parent
        if ~isempty(parents)
            jLL = -inf(1, G) ;
            for j = parents
                X = [dataFix{node} yRHS(:, setdiff(parents, j))];
                jmdl = fitMdl(X, yLHS(:, node)) ;
                jLL(j) = jmdl.LogLikelihood ;                 
            end
            [iLL, iParent] = max(jLL) ; 
            
            % See if new solution increases BIC score. 
            ilogL = logL ;
            ilogL(node) = iLL ;
            nParams = cnt + sum(sum(inter)) - 1 ; 
            iBIC = mdlScore(sum(ilogL), nParams, nObs) ;
        
        
            % See if BIC is lower, if lower approve
            if iBIC < BIC 
                inter(iParent, node) = 0 ; % inter(PARENT, NODE)
                BIC = iBIC ;
                logL = ilogL ;
                %disp(['Gene ', num2str(node), ' removed parent ', num2str(iParent)])
            end
        end
    end
end
end

function score = mdlScore(logL, cnt, nObs)
% computes scores for BIC,EB and VB. Takes scalar values as input. 
% score should be as low as possible to indicate good fit. The bayesian
% score is therefore score = -logEvidence
    para = get_parameters() ;

    if strcmpi(para.robustOpts, 'gaussian') || strcmpi(para.robustOpts, 'constML')
        % [aic,bic] = aicbic(sum(logL), cnt, nObs) ; 
        score = -2*logL + cnt.*(log(nObs) - log(2*pi)) ;

    elseif strcmpi(para.robustOpts, 'EB') || strcmpi(para.robustOpts, 'VB')
        % mdl gives logEvicence as logL with Bayesian Est.
        score = -logL ;
    end
       

end


