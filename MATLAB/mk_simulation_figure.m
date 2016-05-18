function mk_simulation_figure(logData, robustOpts)
% This function makes the simulations of the Drosophila. It can easily be
% changed from ML to Baysian by changing robustOpts and the top captions of
% the figures. 

% Load parameters
para = get_parameters() ;

% Passed arguments
para.robustOpts = robustOpts; % constML, VB

% HardCode parameters of experiment
para.diff_type = 'euler' ; 
para.endo = [1 2 3 4] ;
para.nSim = 1000;

% Set up networks 
interCell{1} = get_network('RPJ') ;
interCell{2} = get_network('ST') ;
interCell{3} = get_network('JAGR') ;
interCell{4} = get_network('PJRG') ;
interCell{5} = get_network('LRN') ; 
[interCell{6}, genes] = get_network('FULL') ;


names = {'RPJ','ST','JAGR','PJRG','LRN','FULL'} ;

% Predim
nCol = length(interCell) + 1 ; % number of networks analyzed + target values
pred = cell(1, nCol-1) ;
sigma =  cell(1, nCol-1) ;
bic = cell(1, nCol-1) ;

% Simulate networks
for c = 1:nCol-1
    [pred{c}, bic{c}] = mk_simulation(logData, interCell{c}, para) ;
    sigma{c} = squeeze(mean(std(data_trf(pred{c}, 'bkw'), 0, 4), 2)) ; % sigma(Gene,Time)
    pred{c} = data_trf(mean(pred{c}, 4), 'bkw') ;
end

% Crate Target data from log-data
target = data_trf(logData, 'bkw') ;

% make figure
figure % initialize something

loc = 35:92 ; % spatial dimension 
tVec = [1 3 5 7 8] ; % which frames to include
nRow = length(tVec) ;
pos = 1:nRow*nCol ;
pos = reshape(pos, nCol, nRow)' ;

% simData LHS
for c = 1:nCol
    if c < nCol
        for r = 1:nRow
            subplot(nRow, nCol, pos(r,c)) ; %
            plot(loc, squeeze(pred{c}(:,:,tVec(r)))') 
            
            % Set dimensions and heading
            set(gca,'ylim',[0 255])
            set(gca,'xlim',[35 92])
            ylabel(['t = ' num2str(tVec(r))])
            if r == 1
                title([names{c},' (', num2str(round(bic{c}, 2)),')'] )
            end
            % Don't add confidence bounds to bayesian methods
            if sum(strcmp(robustOpts,{'VB', 'EB'} )) == 0 
                hold on 
                plot(loc, squeeze(pred{c}(para.endo,:,tVec(r)))' + repmat(1.96*sigma{c}(para.endo, tVec(r)),1 , length(loc))',':r')
                plot(loc, squeeze(pred{c}(para.endo,:,tVec(r)))' - repmat(1.96*sigma{c}(para.endo, tVec(r)),1 , length(loc))',':r')
                hold off
            end
        end
    else
        for r = 1:nRow
            subplot(nRow, nCol, pos(r,c)) ; %
            plot(loc, squeeze(target(:,:,tVec(r)))') 
            set(gca,'ylim',[0 255])
            set(gca,'xlim',[35 92])
            ylabel(['t = ' num2str(tVec(r))])
            if r == 1
                title('TRGT')
            end
        end
    end
end

% Create legend with names
hL = legend(genes,'Location','southoutside','Orientation','horizontal') ;
set(hL,'Position', [0.4 -0.05 0.2 0.2],'Units', 'normalized');    
end


% 
%         subplot(nRow, 3, c2(r)) ; %
%         plot(loc, squeeze(pred2(:,:,tVec(r)))') 
%         set(gca,'ylim',[0 255])
%         set(gca,'xlim',[35 92])
%         ylabel(['t = ' num2str(tVec(r))])
%         if r == 1
%             %title(['Learned BIC:', num2str(bic2)])
%             title(['Learned Ev:', num2str(bic2)])
%         end
% 
%         subplot(nRow, 3, c1(r)) ; %
%         plot(loc, squeeze(realData(:,:,tVec(r)))') 
%         set(gca,'ylim',[0 255])
%         set(gca,'xlim',[35 92])
%         ylabel(['t = ' num2str(tVec(r))])
%         if r == 1
%             title('Target Data')
%         end
% 
