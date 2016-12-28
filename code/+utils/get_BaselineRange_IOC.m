function BaseLineRange= get_BaselineRange_IOC(unit_type)
    
    if nargin <1 
        unit_type = 'samples';
    end
    
    
    WindowOfInterest    = utils.get_WindowOfInterest_IOC('sample');     %  5    60
    BaseLineRange       = 100+[(5*5),(17*5)]; %5kHz Fs!  5    17
        
    if regexpi(unit_type,'sample')
        BaseLineRange       = [(1+BaseLineRange(1)-WindowOfInterest(1)),BaseLineRange(2)-WindowOfInterest(1)];
    elseif regexpi(unit_type,'ms')
        BaseLineRange       = (BaseLineRange-WindowOfInterest(1))./5; %5kHz Fs!        
    else 
        error('getBSL:UnitInvalid','No Valid unit');
    end
    
    
    
end

