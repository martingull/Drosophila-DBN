function data = data_trf(data, direction, para) 
% Applies a transform to data. 
% direction: 'fwd' - transform data, 'bwd' - inverse transformed data

if nargin < 3
    para = get_parameters() ;
end

if strcmpi('log_50', para.transform)

    if strcmpi('fwd', direction)
        data = log((1 + data)/50) ;
        
    elseif strcmpi('bkw', direction)
        data = exp(data) * 50 - 1 ; 
        
    end
    
elseif strcmpi('log256', para.transform)

    if strcmpi('fwd', direction)
        data = log((1 + data)/256) ;
        
    elseif strcmpi('bkw', direction)
        data = exp(data) * 256 - 1 ; 
        
    end
    
elseif strcmpi('none', para.transform)
    
    % return untransformed data
    
else    
    
    disp('Transform in get_paramters() not supported')
    
end

end

