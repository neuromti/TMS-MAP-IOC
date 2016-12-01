function [HS] = get_GroupHS(SUB)
% @params facultatory SUB
% @returns HS, vector of average hotspot
%               if SUB has been used as parameter, presents individual´s HS
% Based on the average across the group  
  % grpstats(squeeze(cat(1,SUB.HS)),SUBID)

  if nargin>0
      HS = single(grpstats(squeeze(cat(1,SUB.HS)),cat(1,SUB.subID)));
  else
    % literals for speed
    % returning the average
    HS = single([ -47.5430  -15.9830   89.0851]);
  end

  
end
