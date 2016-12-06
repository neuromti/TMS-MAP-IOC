% CONSTRUCTION OF THE PERMUTATION ARRAY
% Permutes the factor levels for all subjects present in the true data 
% @returns PERM, a subject-balanced permutation matrix 
function PERM = get_PermMatrix(SUBID,NUM_REP)
    PERM        = int32([]);
    u_sub       = unique(SUBID);
    for rep=1:NUM_REP, 
        perm_set    = NaN(size(SUBID));
        for i_sub=1:length(u_sub)
            sub_GetPut  = find(SUBID==u_sub(i_sub));      
            perm_sample = sub_GetPut(randperm(length(sub_GetPut)),:);
            perm_set(sub_GetPut) = perm_sample;
        end
        %output test all((SUBID-SUBID(perm_set))==0)
        PERM(:,rep) = perm_set;
    end
end