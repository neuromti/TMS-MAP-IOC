function [mapping,sub] = prepare_data(mapping,sub)    

% Takes only values form the initial 5x7 Grid
% Calculates MEP+, centers Latency
mapping.amp                     = single(mapping.amp(1:315,:));
mapping.lat                     = single(mapping.lat(1:315,:));
mapping.xyz                     = single(mapping.xyz(1:315,:));
mapping.mep                     = utils.calculate_MEP(mapping.amp,mapping.lat);
mapping.lat                     = utils.center_latency(mapping.lat);

% Removes artifacted trials, calculates MEP+, centers Latency
remove_flag                     = utils.find_artifactedTrials(mapping);
sub.xyz                         = mapping.xyz(~remove_flag,:);
sub.amp                         = mapping.amp(~remove_flag);    
sub.lat                         = utils.center_latency(mapping.lat(~remove_flag));
sub.mep                         = utils.calculate_MEP(sub.amp,sub.lat);

end
