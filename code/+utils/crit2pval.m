% Calculation of permutation-derived P-Values
% @params Local_Test_Sval, the statistical value from the real experiment
% @params Local_Perm_Sval, a vector of statistical values from the permutated experiments
% @returns P, an estimate of the alpha error
function P = crit2pval(TrueSval,RandomSval)
    TrueMat     = repmat(TrueSval,size(RandomSval,2),1)';
    numSig      = sum(RandomSval>=TrueMat,2);
    divideBy    = (size(RandomSval,2));
    P           = single(numSig./divideBy);
    P           = max(P,1/divideBy)'; 
end