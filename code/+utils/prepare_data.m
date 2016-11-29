function [mapping] = prepare_data(mapping)    

% Takes only values form the initial 5x7 Grid
% Calculates MEP+, centers Latency
mapping.IsArtifacted                = utils.find_artifactedTrials(mapping);
mapping.lat(mapping.IsArtifacted)   = NaN;         
mapping.mep(mapping.IsArtifacted)   = NaN;
mapping.mep                         = utils.calculate_MEP(mapping.amp,mapping.lat);
mapping.lat                         = utils.center_latency(mapping.lat);

end
