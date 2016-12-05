function ProcessingTime = measure_ProcessingTime(starttime)
    if nargin <1 ,
        ProcessingTime  = datetime('now');
    elseif isdatetime(starttime)
        O               = datetime('now')-starttime;
        ProcessingTime  = sprintf('It ran for %.1g years %.1g months %.1g days %.1g hours %.2g minutes %.2g seconds',datevec(O));
    else
        error('Not a valid datetime object');
    end
        

end
