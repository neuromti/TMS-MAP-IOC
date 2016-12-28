function IsArtifacted = find_ArtifactedSubjects(lcl_folder,D)
% @param lcl_folder is the folder where datasets are stored (obligatory)
% @param D is the files in the lcl_folder (facultativ)
% @return IsArtifacted, a logical vector with [1] indicating artifacts
if nargin <2
    D           = dir(lcl_folder);
    D([1,2])    = [];
end

IsArtifacted = true(size(D));
for idx_dataset = 1:length(D)
    if regexp(D(idx_dataset).name,'\w*!!!.mat')%artifacted sets were marked with three exclamation marks
        continue;
    end
    load([lcl_folder,D(idx_dataset).name]);    
    if ~any(mapping.amp>0)           
        continue;
    end           
    IsArtifacted(idx_dataset) = false;
end
