function TestResults = TimeCourse2statistics(DATA,DESIGN,SUBID,STIM)

%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
% -------------------------------------------------------------------------
clear TestResults
WindowOfInterest    = utils.get_WindowOfInterest_IOC('sample');
TimePoints          = (WindowOfInterest(1):WindowOfInterest(2));
if size(DATA,2)~=length(TimePoints), error('Temporal Mismatch between Data and WindowOfInterest'); end

delete(gcp('nocreate'));
parpool(4)
parfor i_timepoint=1:length(TimePoints), 
    TimePointData                           = squeeze(DATA(:,i_timepoint,:));    
    tmp                                     = utils.iocdata2statistics(TimePointData,DESIGN,SUBID,STIM,true,true);
    tmp.Time_ms                             = (TimePoints(i_timepoint)-100)./5;  %transform sample to ms
    TestResults(i_timepoint)                = tmp;       
end
delete(gcp('nocreate'));


