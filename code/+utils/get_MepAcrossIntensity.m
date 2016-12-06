function MepAcrossIntensity = get_MepAcrossIntensity(ioc)

    % Prepare data
    AmplitudeAsMatrix       = reshape(ioc.amp,10,7);        
    LatencyAsMatrix         = reshape(ioc.lat,10,7);
    
    % Find artifacted trials and trials without MEP
    
    IsArtifactAsMatrix  = reshape(utils.find_Artifacted_IOC(ioc),10,7);
    IsMEPAsMatrix       = reshape(ioc.lat==0,10,7);
    
    % Remove artifacted trials and trials without MEP
    LatencyAsMatrix((IsMEPAsMatrix | IsArtifactAsMatrix)) = NaN;
    AmplitudeAsMatrix(IsArtifactAsMatrix) = NaN;
    
    % Threshold MEPs
    % Mills-Nithi
    % MepAsMatrix = (AmplitudeAsMatrix>20) & (LatencyAsMatrix>17) & (LatencyAsMatrix<30);
    
    % Rossini-Rothwell: 
    MepAsMatrix = (AmplitudeAsMatrix>50) & (LatencyAsMatrix>0);
    

    % Calculate Average across intensities
    MepAcrossIntensity = nanmean(MepAsMatrix,1)';

end