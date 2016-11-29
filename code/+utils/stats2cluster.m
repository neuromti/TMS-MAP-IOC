function [ClusVal,GridIdx] = stats2cluster(Hgrid,Sgrid)

[B,GridIdx,N,~]   = bwboundaries(Hgrid,'noholes');
ClusVal     = [];
%ClusSize    = [];
if ~isempty(B),
    for i_clus = 1:N,
        ClusVal     = [ClusVal,sum(Sgrid(GridIdx==i_clus))];
    %    ClusSize    = [ClusSize,sum(sum(GridIdx==i_clus))];
    end

    ClusVal     = sort(ClusVal,'descend');
  %  ClusSize    = sort(ClusSize,'descend');

else
    ClusVal     = [];
   % ClusSize    = [];
end
end