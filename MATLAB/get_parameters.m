function para =  get_parameters()
% Default settings for simulation.
% These are typically changed within each experiment.


% Fixed Effects
para.autoreg = 1 ; % 1 - first order autoreg, 0 - no autoreg
para.diff_type = 'euler' ; % see mk_euler for options,

% Endogenous Columns
para.endo = [1 2 3 4] ; % columns of endo genes, exs [1:4]

% Transform Choice
para.transform = 'log256' ; 

% Set search values.
para.nIter = 200 ; % Evaluations of causal structure
para.nRestart = 4 ; % 
para.robustOpts = 'VB' ; % see fitMdl, node assumptions. 

% Set simulations parameters
%para.mkSim = 'on' ;
para.nSim = 1000 ; % Number of times to simulate system

% Set movie parameters
para.movie = 'off' ; % Make move: 'on','off'
para.nFrames = 180 ; % Total frames 
para.frameRate = 30 ; % Frames per second


end