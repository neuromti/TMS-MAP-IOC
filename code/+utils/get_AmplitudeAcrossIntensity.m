function AmplitudeAcrossIntensity = get_AmplitudeAcrossIntensity(ioc)

    % Prepare data
    AmplitudeAsMatrix     = reshape(ioc.amp,10,7);
    
    % Find artifacted trials and trials without MEP
    
    IsArtifactAsMatrix      = reshape(utils.find_Artifacted_IOC(ioc),10,7);
    
    % Remove artifacted trials and trials without MEP
    AmplitudeAsMatrix(IsArtifactAsMatrix) = NaN;
    
    % Calculate Average across intensities
    AmplitudeAcrossIntensity = nanmean(AmplitudeAsMatrix,1)';

end