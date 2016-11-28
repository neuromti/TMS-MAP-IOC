function progressBar(current_rep)
TrailArray      = {'-','\','|','/','-','\','|','/'};
MOD_CONSTANT    = 10;


if ischar(current_rep)
    % First
    if strcmpi(current_rep,'0'),
        fprintf('['),
        fprintf('%s',TrailArray{1});
    %  Lasz
    elseif strcmpi(current_rep,'1'),
        fprintf('\b%s','.')   
        fprintf('%s\n',']')       
    else
        warning('Error, Use "0" and "1" to mark Start and End');
    end
    
    return
end    

% Running
if ~(mod(current_rep,MOD_CONSTANT)),          
    tl = mod(int32(current_rep/MOD_CONSTANT),length(TrailArray))+1;        
    if tl==1,
        fprintf('\b%s','.')            
        fprintf('%s',TrailArray{tl});
    else
        fprintf('\b%s',TrailArray{tl});
    end  
else
    return;
end

end