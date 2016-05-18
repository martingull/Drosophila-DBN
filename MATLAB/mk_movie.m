function mk_movie(data, movie_name)
    % Takes a 3 dimentional vector (gene x space x time)
    % The vector is interpolated in time (makes movie smooth)
    % Movie is then saved as MP4 format
    
    para = get_parameters() ;
    if strcmpi(para.movie, 'on')
    
    % Number of Frames
    nFrames = para.nFrames ;
    frameRate = para.frameRate ;
    
    % Interpolate the input
    t = 1:size(data,3) ;
    tq = linspace(1, size(data,3), nFrames) ; 
   
    newData = zeros(size(data,1), size(data,2), nFrames) ;
    for g = 1:size(data,1)
        for i = 1:size(data,2)
            newData(g,i,:) = spline(t, squeeze(data(g,i,:))', tq) ;
        end
    end
    
    % Setting Up Movie writer 
    writerObj = VideoWriter(movie_name, 'MPEG-4') ;
    writerObj.FrameRate = frameRate ;
    open(writerObj) ;
    
    axis tight
    %legend(para.names) % don't know set command for legend. 
    set(gca,'nextplot','replacechildren');
    set(gcf,'Renderer','zbuffer');
    set(gca,'ylim',[0 255] )
    
    % Write each frame of the Interpolated data
    for t = 1:size(newData,3)
        plot(squeeze(newData(:,:,t))') 
        %
        frame = getframe ; 
        writeVideo(writerObj, frame)
    end

    close(writerObj)
    else
        % Don't make a movie
    end
 
end