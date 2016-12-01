function [ANT] = get_GroupAnt(SUB)
% @params facultatory SUB
% @returns ANT, vector of average hotspot
%               if SUB has been used as parameter, presents individual´s ANTs
% Based on the average across the group  
  % grpstats(squeeze(cat(1,SUB.ANT)),SUBID)

  if nargin>0
    K           = cat(1,SUB.subID);
    HS          = grpstats(squeeze(cat(1,SUB.HS)),K,{'mean'});
    ANT         = grpstats(squeeze(cat(1,SUB.ANT)),K,{'mean'});
    ANT         = ANT-HS+(repmat(utils.get_DesignGridOrigin+[0,-10,0],size(HS,1),1));  
  else
    % literals for speed
    % returning the average
    ANT = single([-32.3619    3.4596   66.8136]);    
  end

  
end
