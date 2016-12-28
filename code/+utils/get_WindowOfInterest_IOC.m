function WindowOfInterest = get_WindowOfInterest_IOC(unit_type)
    if nargin <1
        unit_type = 'samples';
    end
    
    if regexpi(unit_type,'sample')
        WindowOfInterest =  100+[(5*5),(60*5)];
    elseif regexpi(unit_type,'ms')
        WindowOfInterest =  ([(5*5),(60*5)])./5;
    else 
        error('WOI:UnitValidity','No Valid unit');
    end
end