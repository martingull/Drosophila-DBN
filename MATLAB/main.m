clear 

% Global Parameters, set in get_parameters()
para = get_parameters() ;


% Load data, transform data
load data_Nankang.mat % Nankang Data, (Gene,Instance,Time)
data = data(:,:,2:end) ; % First time frame from C13
logData = data_trf(data, 'fwd') ;


% Searching for causal structure.
inter = search_dbn(logData) ; 
%inter = get_network('LRN') ; % Precomputed

% Make Simulation Proposed Network (Makes Picture)
mk_simulation_figure(logData, 'constML')

% Make matrix of RMSE
RMSE = performance_analysis(data, 'constML') ;  

% Make matrix of Regulatory weights
regWgt = mk_regulation_matrix(data, 'constML') ;

% Experiment: Make list of probable networks, 
%mk_experiment_metropolis_edges(logData)  ;

% Make movie
%mk_movie(simData, 'sim_Movie')

