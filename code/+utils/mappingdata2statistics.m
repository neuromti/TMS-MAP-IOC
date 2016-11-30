function [TestResults,ClusterResults] = mappingdata2statistics(DATA,SUBID,DESIGN)
%% -------------------------------------------------------------------------
% DEFINITION OF PARAMETERS
% -------------------------------------------------------------------------

ALPHA_ERROR         = 0.05;
DESIGN_MATRIX       = cat(2,DESIGN,SUBID);
NUM_REP             = 1000;
disp(['Alpha Error: ',num2str(ALPHA_ERROR)])
disp(['Number of Repetitions ',num2str(NUM_REP)])
%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
% -------------------------------------------------------------------------

for i_pos=1:size(DATA,1),     
    PositionData                    = (DATA(i_pos,:))';
    [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(PositionData,DESIGN_MATRIX);       
    Test_Pval(i_pos,:)              = Pval;        
    Test_Sval(i_pos,:)              = StatVal;
    Test_Coeffs(i_pos,:)            = Coeffs;
    Test_MargMeans(i_pos,:)         = MargMeans;
end

for k=1:size(Test_Sval,2),    
    Hgrid                       = utils.vec2mesh((Test_Pval(:,k)<ALPHA_ERROR));
    Sgrid                       = utils.vec2mesh(Test_Sval(:,k));    
    [ClusVal,GridIdx,ClusSize]  = utils.stats2cluster(Hgrid,Sgrid);    
    Test_ClusVal{k}             = ClusVal;
    Test_ClusSize{k}            = ClusSize;
    Test_ClusIdx{k}             = GridIdx;
end
       
TestResults              = struct('Pval',Test_Pval,'Sval',Test_Sval,'Coeffs',Test_Coeffs,'MargMeans',Test_MargMeans);
ClusterResults           = struct('Threshold',ALPHA_ERROR,'MemberShip',Test_ClusIdx,'Size',Test_ClusSize,'Sval',Test_ClusVal);

%% ------------------------------------------------------------------------
% DISTRIBUTION FREE TEST
% -------------------------------------------------------------------------
% Check whether they should be performed
if NUM_REP>0 
    disp('Starting permutation and bootstrap analysis!')
else
    disp('You chose not to run permutation and bootstrap analysis!')
    return;
end

% PERMUTATION ANALYSIS 
PERM    = get_PermMatrix(SUBID,NUM_REP);
for i_pos=1:size(DATA,1),        
    utils.progressBar(['Permutation GridPoint ',num2str(i_pos),' [.']);    
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        PermutationData                 = (DATA(i_pos,PERM(:,rep)));          
        [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(PermutationData,DESIGN_MATRIX);          
        Perm_Pval(i_pos,:,rep)          = Pval;        
        Perm_Sval(i_pos,:,rep)          = StatVal;
        Perm_Coeffs(i_pos,:,rep)        = Coeffs;
        Perm_MargMeans(i_pos,:,rep)     = MargMeans;
    end    
    testCriterion   = abs(squeeze(Test_Coeffs(i_pos,:,:)));
    permCriterion   = abs(squeeze(Perm_Coeffs(i_pos,:,:))); 
    Perm_ContrastPval(i_pos,:)          = perm2pval(testCriterion,permCriterion);
    utils.progressBar('1');
end
TestResults.PermutationPval = Perm_ContrastPval;


% BOOTSTRAP ANALYSIS 
BOT     = get_BootMatrix(DESIGN,NUM_REP);
for i_pos=1:size(DATA,1),        
    utils.progressBar(['Bootstrap GridPoint ',num2str(i_pos),' [.']);    
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        BootData                        = (DATA(i_pos,BOT(:,rep)));          
        [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(BootData,DESIGN_MATRIX);          
        Boot_Pval(:,rep)          = Pval;        
        Boot_Sval(:,rep)          = StatVal;
        Boot_Coeffs(:,rep)        = Coeffs;
        Boot_MargMeans(:,rep)     = MargMeans;
    end
    utils.progressBar('1');   
    Boot_Power(i_pos,:)     = nanmean(Boot_Pval<ALPHA_ERROR,2);
    Boot_CoeffsCI(i_pos,:,:)= cat(2,quantile(Boot_Coeffs,ALPHA_ERROR./2,2),quantile(Boot_Coeffs,(1-ALPHA_ERROR./2),2));      
end
TestResults.Power       = single(Boot_Power);
TestResults.CoeffsCI    = single(Boot_CoeffsCI);


%% ------------------------------------------------------------------------
% PERMUTATION BASED CLUSTER ANALYSIS 
% -------------------------------------------------------------------------

for k=1:length(ClusterResults)
    TestClusterNum                      = length(ClusterResults(k).Sval);    
    
    if TestClusterNum==0,    
        ClusterResults(k).Perm_Pval = [];
        continue; 
    end            
    
    if TestClusterNum>0  
        
        [AbsTestClusterVal,WhichTestClus]   = sort(abs(ClusterResults(k).Sval),'descend');   
        [~,PutBackInOrder]                  = sort(WhichTestClus,'ascend');
        SigCounts                           = true(NUM_REP,TestClusterNum);
        
        for rep=1:NUM_REP,        
            Hgrid               = utils.vec2mesh(Perm_Pval(:,k,rep)<ALPHA_ERROR);
            Sgrid               = utils.vec2mesh(Perm_Sval(:,k,rep));
            ClusVal             = utils.stats2cluster(Hgrid,Sgrid);

            if isempty(ClusVal),
                IsSignificantbyChance   = false(1,TestClusterNum);
                SigCounts(rep,:)        = IsSignificantbyChance;  
                continue; 
            end
            
            AbsPermClusterVal   = sort(abs(ClusVal));                        
            
            for PickTestCluster=1:TestClusterNum
                IsSignificantbyChance(PickTestCluster) = AbsPermClusterVal(1)>AbsTestClusterVal(PickTestCluster);
            end
            SigCounts(rep,:)    = IsSignificantbyChance(PutBackInOrder);              
        end               
        tmpPval         = sum(SigCounts)./(NUM_REP);
        ClusterPermPval = max(tmpPval,1/NUM_REP);   
    end    
    
end
   
end
%% LOCAL FUNCTIONS
% ANALYSIS FUNCTION 
% Calculates the  statistical Values
% @params DATA, a vector of data to perform the analysis
% @params DESIGN_MATRIX, the desoing matrix of the experiment
% @return PERM, a subject-balanced permutation matrix 

function [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(DATA,DESIGN_MATRIX)   
    % DEF_MODEL will returning Interaction 
    % Explained as Crosstable:
    %     90  45              90 45
    % BI   +  -     oder  BI  -   +
    % MO   -  +           MO  +   -
    % We show the sign of BI 90° and MO 45°   
    DEF_MODEL               = [1 0 0;0 1 0;1 1 0;0 0 1];        
    % Picking Values will be for 'Waveform', 'Orientation', 'Interaction')
    PickVal                 = 2:4;
    % Picking Coeffs will be for ('BI=1', 'LM=1', 'BI=1 * LM=1')
    PickCoeffs              = [3,5,9];
    statistical_test        = @(x,y)anovan(x,y,'random',3,'display','off','varnames',{'BI','LM','SUB'},'model',DEF_MODEL);    
    [~,tab,stats]           = statistical_test(DATA,DESIGN_MATRIX);
    warning('off','stats:multcompare:IgnoringRandomEffects')
    [~,m,~,~]               = multcompare(stats,'dim',[1,2],'display','off');
    warning('on','stats:multcompare:IgnoringRandomEffects')
    Pval                    = single([tab{PickVal,7}]);
    StatVal                 = single([tab{PickVal,6}]);    
    MargMeans               = single(m(:,1));
    Coeffs                  = single(stats.coeffs(PickCoeffs));
end

% CONSTRUCTION OF THE PERMUTATION ARRAY
% Permutes the factor levels for all subjects present in the true data 
% @returns PERM, a subject-balanced permutation matrix 
function PERM = get_PermMatrix(SUBID,NUM_REP)
    PERM        = int8([]);
    u_sub       = unique(SUBID);
    for rep=1:NUM_REP, 
        perm_set    = [];
        for i_sub=1:length(u_sub)
            sub_GetPut  = find(SUBID==u_sub(i_sub));        
            perm_set    = cat(1,perm_set,sub_GetPut(randperm(length(sub_GetPut)),:));
        end
        PERM(:,rep) = perm_set;
    end
end

% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Draws randomly from subjects with replacement for all factor levels present in the true data 
% @returns BOT, a factor-balanced bootstrap matrix
function BOT = get_BootMatrix(DESIGN,NUM_REP)
    BOT         = int8([]);
    dec_design  = bin2dec(num2str(DESIGN))+1;
    u_dsg       = unique(dec_design);
    for rep=1:NUM_REP, 
        bot_set    = [];
        for i_dsg=1:length(u_dsg)
            dsg_GetPut  = dec_design==u_dsg(i_dsg);
            bot_set     = cat(1,bot_set,datasample(find(dsg_GetPut),sum(dsg_GetPut)));
        end
        BOT(:,rep) = bot_set;
    end
end


% Calculation of permutation-derived P-Values
% @params Local_Test_Sval, the statistical value from the real experiment
% @params Local_Perm_Sval, a vector of statistical values from the permutated experiments
% @returns P, an estimate of the alpha error
function P = perm2pval(TrueSval,RandomSval)
    TrueMat     = repmat(TrueSval,size(RandomSval,2),1);
    numSig      = sum(RandomSval'>=TrueMat);
    divideBy    = (size(RandomSval,2));
    P           = single(numSig./divideBy);
    P           = max(P,1/divideBy); 
end