function progressBar(current_rep,MOD_CONSTANT)

TrailArray      = {'-','\','|','/','-','\','|','/'};

if nargin<2,  MOD_CONSTANT    = 10; end


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
        fprintf('%s',current_rep)   
        %warning('Error, Use "0" and "1" to mark Start and End');
    end
    
    return
end    

% Running Index
if isnumeric(current_rep)
    if ~(mod(current_rep,MOD_CONSTANT)),          
        tl = mod(int32(current_rep/MOD_CONSTANT),length(TrailArray)+1)+1;        
        pause(1)
        if tl>length(TrailArray),
            fprintf('\b%s','.')                        
        elseif tl==1,
            fprintf('%s',TrailArray{tl});
        else
            fprintf('\b%s',TrailArray{tl});
        end  
    else
        return;
    end
end

end
