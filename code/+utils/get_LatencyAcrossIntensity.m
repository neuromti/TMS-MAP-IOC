function LatencyAcrossIntensity = get_LatencyAcrossIntensity(ioc)

    % Prepare data
    LatencyAsMatrix     = reshape(ioc.lat,10,7);
    
    % Find artifacted trials and trials without MEP
    
    IsArtifactAsMatrix  = reshape(utils.find_Artifacted_IOC(ioc),10,7);
    IsMEPAsMatrix       = reshape(ioc.lat==0,10,7);
    
    % Remove artifacted trials and trials without MEP
    LatencyAsMatrix((IsMEPAsMatrix | IsArtifactAsMatrix)) = NaN;
    
    % Calculate Average across intensities
    LatencyAcrossIntensity = nanmean(LatencyAsMatrix,1)';

end