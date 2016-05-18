function [inter, names] = get_network(net)
% returns intercellular adjacency matrix.

% Assumption: Hkb was included to regulate anterior section of Hb. Adjacency
% from Hkb -> Hb is assumed for all networks. 

if strcmp(net, 'FULL') % Fully connected network. 
    inter = [ones(8, 4) zeros(8, 4)] ;
    inter(logical(eye(size(inter)))) = 0 ;
    
    
elseif strcmp(net, 'RPJ') % Rivera-Pomar Jackle
    inter = [0 1 0 0 1 0 0 1 ; 1 0 1 1 1 0 1 0 ; 0 1 0 1 1 1 1 0 ; 1 1 1 0 1 1 1 0]' ;
    inter = [inter zeros(8, 4) ] ;
    %inter = [0 1 1 1 1 1 1 1 ; 1 0 1 1 1 0 1 0 ; 0 1 0 1 1 1 1 0 ; 1 1 1 0
    %1 1 1 0 ]' ; % Full Hb regularisation.
    
    
elseif strcmp(net, 'ST') % Sanchez Thieffy
    inter = [0 1 0 0 1 0 0 1 ; 1 0 1 1 1 0 0 0 ; 1 1 0 0 1 1 0 0 ; 1 0 1 0 1 1 0 0]' ;
    inter = [inter zeros(8, 4) ] ;
    
    
elseif strcmp(net, 'JAGR') % Jaeger et al.
    inter = [0 0 1 1 1 1 0 1 ; 1 0 1 1 1 1 1 0 ; 1 1 0 0 1 1 1 0 ; 1 1 1 0 1 1 1 0]' ;
    inter = [inter zeros(8, 4) ] ;
    
    
elseif strcmp(net, 'PJRG') % Perkins Jaeger Reinitz Glass
    inter = [0 1 0 1 1 0 1 1 ; 1 0 1 1 1 0 1 0 ; 1 1 0 0 1 1 1 0 ; 1 0 1 0 1 1 1 0]' ;
    inter = [inter zeros(8, 4) ] ;
    
elseif strcmp(net, 'LRN') % The Learned Network
    inter = [0 1 1 1 0 1 1 1 ; 0 0 1 1 0 1 1 0 ; 1 0 0 0 1 0 1 0 ; 1 1 0 0 1 0 1 0]'  ;
    inter = [inter zeros(8, 4)] ;
    
end


names = {'Hb', 'Kr', 'Gt', 'Kni', 'Bcd', 'Cad', 'Tll', 'Hkb'}' ;

end