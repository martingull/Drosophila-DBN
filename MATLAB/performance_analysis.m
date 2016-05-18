function results = performance_analysis(data, robustOpts) 
% This function creates the RMSE table of the thesis. 
% 

% get parameters and hardcode experiment
para = get_parameters() ;

% Hardcode experiment specifics
para.robustOpts = robustOpts ;
para.nSim = 5000 ;
logData = data_trf(data, 'fwd') ;


% Set Nets
interCell{1} = get_network('RPJ') ;
interCell{2} = get_network('ST') ;
interCell{3} = get_network('JAGR') ;
interCell{4} = get_network('PJRG') ;
interCell{5} = get_network('LRN') ; 
interCell{6} = get_network('FULL') ;


results = cell(1, 2) ;
for p = 1:size(interCell, 2) % For each causal netowrk
        
        % Simulate Data using 
        [predData, score]  = mk_simulation(logData, interCell{p}, para) ;
    
        % Transform data to original scale
        predData = data_trf(predData, 'bkw') ;
        simData = mean(predData ,4) ;
    
        for e = 1:4 ;
            results{1}(p,e) = squeeze(sqrt(mean(mean((data(e,:,:)-simData(e,:,:)).^2)))) ;
            results{2}(p,e) = squeeze(sqrt(mean((data(e,:,end)-simData(e,:,end)).^2))) ;
        end
        results{1}(p,5) = score ;
        
end


end