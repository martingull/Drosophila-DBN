function mk_figure(realData, predData, inter, para)
% Todo As number of figures increases make them into functions here.

if nargin<4
    para = get_parameters() ;
    endo = para.endo ;
end
loc = 36:91 ; % Check !

simData = mean(predData, 4) ;

% make confidence bounds plot. 
if sum(strcmpi('conf', para.simFigure)) > 0

sigma = squeeze(mean(std(predData, 0, 4), 2)) ; % sigma(Gene,Time) 

figure % initialize something

tVec = [3 5 7 8] ; % which frames to include
h = length(tVec) ;


% simData LHS
for i = 1:h
    subplot(h, 1, i) ; %
    plot(loc, squeeze(simData(:,:,tVec(i)))')  
    hold on
    %sigma = std(sqeeze(simData(:,:,tVec(i)))') ;
    plot(loc, squeeze(simData(endo,:,tVec(i)))' + repmat(1.96*sigma(endo, tVec(i)),1 , length(loc))',':r')
    plot(loc, squeeze(simData(endo,:,tVec(i)))' - repmat(1.96*sigma(endo, tVec(i)),1 , length(loc))',':r')
    set(gca,'ylim',[0 255])
    set(gca,'xlim',[35 92])
    ylabel(['t = ' num2str(tVec(i))])
    if i == 1
        title('Simulated Confidence Bounds')
    end
    hold off
end
% Create legend with names
hL = legend(para.names,'Location','southoutside','Orientation','horizontal') ;
set(hL,'Position', [0.4 -0.05 0.2 0.2],'Units', 'normalized');        
end




if sum(strcmpi('sim', para.simFigure)) > 0

figure % initialize something

tVec = [1 3 5 7 8] ; % which frames to include
h = length(tVec) ;
left = 1:2:10 ;
right = 2:2:10 ;


% simData LHS
for i = 1:h
    subplot(h, 2, right(i)) ; %
    plot(loc, squeeze(simData(:,:,tVec(i)))') 
    set(gca,'ylim',[0 255])
    set(gca,'xlim',[35 92])
    ylabel(['t = ' num2str(tVec(i))])
    if i == 1
        title('Simulated Data')
    end
    
    subplot(h, 2, left(i)) ; %
    plot(loc, squeeze(realData(:,:,tVec(i)))') 
    set(gca,'ylim',[0 255])
    set(gca,'xlim',[35 92])
    ylabel(['t = ' num2str(tVec(i))])
    if i == 1
        title('Target Data')
    end

end
% Create legend with names
hL = legend(para.names,'Location','southoutside','Orientation','horizontal') ;
set(hL,'Position', [0.4 -0.05 0.2 0.2],'Units', 'normalized');    
end

% Make adjacency Grap
if strcmpi(para.simFigure, 'adjGraph')
    g = drawNetwork('-adjMat', inter, '-nodeLabels', para.names, '-splitLabels', 1, '-layout', Circlelayout);
end


end