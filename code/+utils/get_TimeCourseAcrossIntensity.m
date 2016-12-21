function TimeCourseAcrossIntensity = get_TimeCourseAcrossIntensity(ioc)

    % 5ms after TMS pulse
    % baseline is from 5 to 17ms after TMS pulse
    % Remember 5KHz sampling rate 20ms -> 100 samples
    WindowOfInterest    = utils.get_WindowOfInterest_IOC('sample');
    BaseLineRange       = utils.get_BaselineRange_IOC('sample');
   
    % Prepare Data
    RawWaveform         = reshape(ioc.raw,7*10,[]);
    RawWaveform         = RawWaveform(:,WindowOfInterest(1):WindowOfInterest(2));
    
    % Baseline relative to 20 ms prior to TMS pulse    
    BaselinedWaveform   = utils.baseline(detrend(RawWaveform')',BaseLineRange,1);
   
    % Find artifacted trials and trials without MEP
    IsArtifacted            = utils.find_Artifacted_IOC(ioc);
    CleanedWaveform         = BaselinedWaveform;
    CleanedWaveform(IsArtifacted,:) = NaN;
     
    % Calculate Average across intensities
    WaveformAsMatrix          = reshape(CleanedWaveform,7,10,[]);
    TimeCourseAcrossIntensity = squeeze(nanmean(WaveformAsMatrix,2));
    
    % Smoothen / Gaussian Filter 
    filt_flag  = true;
    if filt_flag == true
        %disp('I smoothen the time course!')
        GaussWinSize    = [2,10];
        filt            = gausswin(GaussWinSize(1))*gausswin(GaussWinSize(2))';
        filt            = filt./sum(sum(filt));
        tmp             = padarray(padarray(TimeCourseAcrossIntensity,[0,GaussWinSize(2)./2],'replicate'),[GaussWinSize(1)./2,0],'replicate');
        output          = convn(tmp,filt,'same');
        output          = output(1+(GaussWinSize(1)/2):end-GaussWinSize(1)/2,1+(GaussWinSize(2)/2):end-GaussWinSize(2)/2);
        TimeCourseAcrossIntensity = output;
    end
%     close all
%     plot(WaveformAcrossIntensity')
    

end