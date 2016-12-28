function [ClusVal,GridIdx,ClusSize] = stats2cluster(Hgrid,Sgrid)

posH                = Hgrid&(Sgrid>0);
negH                = Hgrid&(Sgrid<0);   
[PosClusVal,PosGridIdx,PosClusSize] = calc_cluster(posH,Sgrid);
[NegClusVal,NegGridIdx,NegClusSize] = calc_cluster(negH,Sgrid);

ClusVal     = [PosClusVal,NegClusVal];
ClusSize    = [PosClusSize,NegClusSize];
GridIdx     = cat(3,PosGridIdx,NegGridIdx);

end

function [ClusVal,GridIdx,ClusSize] = calc_cluster(H,S)
    [B,GridIdx,N,~]   = bwboundaries(H,'noholes');
    GridIdx     = int8(GridIdx);
    ClusVal     = single([]);
    ClusSize    = single([]);
    if ~isempty(B)
        for i_clus = 1:N
            ClusVal     = [ClusVal,sum(S(GridIdx==i_clus))];
            ClusSize    = [ClusSize,sum(sum(GridIdx==i_clus))];
        end
    else
        ClusVal     = [];
        ClusSize    = [];
    end

end