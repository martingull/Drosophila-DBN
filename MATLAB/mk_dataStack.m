function [dataFix, yRHS, yLHS] = mk_dataStack(data, para)
% Takes data(gene,instance,time) returns transformed data.

[G, N, T] = size(data) ;
%nObs = (N-2) * (T-1) ;
nObs = N * (T-1) ;


% Allow for custom parameters. 
if nargin < 2
    para = get_parameters() ;
end

% Autoreg Series
yRHS = nan(nObs, G) ; yLHS = nan(nObs, G) ;
for g = 1:G
    if para.autoreg == 1
        %yRHS(:,g) = reshape(data(g,2:end-1,1:end-1), [nObs, 1]) ;
        temp = data(g,:,1:end-1) ;
        yRHS(:,g) = temp(:) ;
        %yLHS(:,g) = reshape(data(g,2:end-1,2:end), [nObs, 1]) ;
        temp = data(g,:,2:end) ;
        yLHS(:,g) = temp(:) ;
    else
	yRHS(g) = [] ;
	yLHS(g) = [] ;
    end
end

% % Squared Autoreg 
% yRHScube = cell(1, G) ;
% for g = 1:G
% 	if para.autoregCube == 1
%         %yRHScube{g} = reshape(data(g,2:end-1,1:end-1).^2, [nObs, 1]) ;
%         temp = data(g,:,1:end-1).^2 ;
%         yRHScube{g} = temp(:) ;
%     else
%         yRHScube{g} = [] ; 
% 	end
% end



% Diffusion Series
diff = cell(1,G) ;
for g = 1:G
    gDiff = mk_euler(squeeze(data(g,:,1:T-1)), para) ;
    if ~isempty(gDiff)
        %diff{g} = reshape(gDiff, [nObs, 1]) ;
        diff{g} = gDiff(:) ;
    else
        diff{g} = [] ;
    end
end

% Concatenate Data
dataFix = cell(1, G) ;
for g = 1:G 
    %dataFix{g} = [ yRHS(:,g) yRHScube{g} diff{g} ] ;
    dataFix{g} = [ yRHS(:,g) diff{g} ] ;
end 

end
