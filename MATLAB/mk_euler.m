function diff = mk_euler(data, para)
% Takes data(space,time) and returns diffusion(space,time).

if nargin < 2
    para = get_parameters() ; % Global parameters 
end

if strcmpi(para.diff_type, 'direct') % direct diffusion
    for t=1:size(data, 2)
        diff = [data(3:end, t) data(1:end-2, t)];
    end
    
elseif strcmpi(para.diff_type, 'euler') % central euler
    diff = nan(size(data)) ;
    diff(1, :) = -data(1, :) + data(2, :) ;
    diff(end, :) = -data(end-1, :)+ data(end, :) ;
    diff(2:end-1,:) = data(1:end-2, :) -2*data(2:end-1, :) + data(3:end, :) ;
    
elseif strcmpi(para.diff_type, 'avg') % average / constrained
    diff = (data(1:end-2, :) + data(3:end, :)).*0.5 ;
    
elseif strcmpi(para.diff_type, 'nghood') % neighbourhood
    
    diff = [data(1, :) ; (data(1:end-2, :) + data(3:end, :)).*0.5 ; data(end, :)] ;
    
elseif strcmpi(para.diff_type, 'off') % no diff
    diff = [] ;
    
elseif strcmpi(para.diff_type, 'perp') % orthogonal 
    diff = [] ;
    for t = size(data,2)
        % X ~ feedback, y ~ diffusion -> eps ~ diffusion without feedback. 
        mdl = fitlm([data(1:end-2, :) data(3:end, :)], data(2:end-1, :)) ;
        diff = [diff mdl.Residuals.Raw] ;
    end
    
elseif strcmpi(para.diff_type, 'vario') % variogram
    diff = 0.5 * (data(1:end-2, :) - data(3:end, :)).^2 ;
    
elseif strcmpi(para.diff_type, 'kernel') % Matern Kernel
    diff = exp(-((data(1:end-2, :) -2*data(2:end-1, :) + data(3:end, :))/0.1).^2) ;
    
elseif strcmpi(para.diff_type, 'kring') % Markov Kringing
    diff = sqrt((data(1:end-2, :)-data(2:end-1, :)).^2 + (data(2:end-1, :)-data(3:end, :)).^2) ;
    
else
    disp('Warning: mk_euler passed incorrect spec.')
    
end

end