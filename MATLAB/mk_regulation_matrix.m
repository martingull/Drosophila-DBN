function outData = mk_regulation_matrix(data, robustOpts) 

% Get Parameters
para = get_parameters() ;


% Hardcode Assumptions
para.robustOpts = robustOpts ;
logData = data_trf(data, 'fwd') ;
[dataFix, yRHS, yLHS] = mk_dataStack(logData, para) ;



% Precomputed Networks
interCell{1} = get_network('RPJ') ;
interCell{2} = get_network('ST') ;
interCell{3} = get_network('JAGR') ;
interCell{4} = get_network('PJRG') ;
interCell{5} = get_network('LRN') ; % Learned Network
interCell{6} = get_network('FULL') ;

nNet = size(interCell, 2) ;
G = 8 ;
E = 4 ;



outData = zeros(nNet, E*G) ;
for i = 1:nNet % For each causal structure
    inter = interCell{i} ;
    
    for j = 1:E % For each endo variable
            parents = find(inter(:, j))' ;
            X = [dataFix{j} yRHS(:, parents)] ;
            
            % Fit the model
            mdl = fitMdl(X, yLHS(:, j), para) ; 
            
            % Returns parameters
            if strcmpi(para.robustOpts, 'gaussian')
                outData(i,(j-1)*8 + parents) = mdl.Coefficients.Estimate(4:end) ; % gaussian 
            elseif strcmpi(para.robustOpts, 'constML')
                outData(i,(j-1)*8 + parents) = mdl.theta(4:end) ; % constrained
            elseif strcmpi(para.robustOpts, 'VB')
                outData(i,(j-1)*8 + parents) = mdl.wN(4:end) ; % VB
            else
                disp('Add more options to mk_regulation_matrix')
            end
    end
end
