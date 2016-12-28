function ALPHA_ERROR  = get_ALPHAERROR(NumTest)
    if nargin == 0
        ALPHA_ERROR         = 0.05;
    elseif nargin == 1 && isnumeric(NumTest)
        ALPHA_ERROR         = 0.05/NumTest;
    else
        error('getALPHA:UnspecifiedTestNumber','No correctly specified number of tests');
    end
    
end