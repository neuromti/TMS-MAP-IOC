function [TestResults,ClusterResults] = iocdata2statistics(DATA,DESIGN,SUBID,STIM,PERM_flag,BOT_flag)
%% -------------------------------------------------------------------------
% DEFINITION OF PARAMETERS
% -------------------------------------------------------------------------
ALPHA_ERROR         = utils.get_ALPHAERROR();
DESIGN_MATRIX       = cat(2,DESIGN,SUBID);
NUM_REP             = 1000;
NUM_STIM            = 7;
disp(['Alpha Error: ',num2str(ALPHA_ERROR)])
disp(['Number of Repetitions ',num2str(NUM_REP)])

% Flags for structured code testing 
if nargin <5,
    PERM_flag   = true;
    BOT_flag   = true;
end
if nargout >1 && PERM_flag
    CLUS_flag   = true;
else
    CLUS_flag   = false;
end

%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
% -------------------------------------------------------------------------

for i_stim=1:NUM_STIM
    StimData                        = DATA(STIM==i_stim,:);
    StimDesign                      = DESIGN_MATRIX(STIM==i_stim,:);
    [Pval,StatVal,Coeffs,MargMeans] = get_StatisticalValues(StimData,StimDesign);       
    Test(i_stim).Pval               = Pval;        
    Test(i_stim).Sval               = StatVal;
    Test(i_stim).Coeffs             = Coeffs';    
    Test(i_stim).MargMeans          = MargMeans;    
end

Test_Pval   = cat(1,Test.Pval);
Test_Sval   = cat(1,Test.Sval);
Test_Coeffs = cat(1,Test.Coeffs);
Test_MargMeans = cat(2,Test.MargMeans)';
TestResults              = struct('Pval',Test_Pval,'Sval',Test_Sval,'Coeffs',Test_Coeffs,'MargMeans',Test_MargMeans);

if CLUS_flag,
for k=1:size(Test_Sval,2),    
    Hgrid                       = Test_Pval(:,k)<ALPHA_ERROR;
    Sgrid                       = Test_Sval(:,k);    
    [TestClusVal,TestGridIdx,TestClusSize]  = utils.stats2cluster(Hgrid,Sgrid);    
    Test_ClusVal{k}             = TestClusVal;
    Test_ClusSize{k}            = TestClusSize;
    Test_ClusIdx{k}             = TestGridIdx;
end
ClusterResults           = struct('Threshold',ALPHA_ERROR,'MemberShip',Test_ClusIdx,'Size',Test_ClusSize,'Sval',Test_ClusVal);
end
%% ------------------------------------------------------------------------
% Check whether they should be performed
if (~PERM_flag && ~ BOT_flag)
    disp('You chose not to run permutation and bootstrap analysis!')
    return;
end

% DISTRIBUTION FREE TEST
% -------------------------------------------------------------------------
if PERM_flag
    disp('Starting permutation analysis!')
% PERMUTATION ANALYSIS 
Perm = struct();
for i_stim=1:NUM_STIM      
    StimData        = DATA(STIM==i_stim,:);
    StimDesign      = DESIGN_MATRIX(STIM==i_stim,:);
    PermSUBID       = SUBID(STIM==i_stim);
    PERM            = get_PermMatrix(PermSUBID,NUM_REP);
    
    utils.progressBar(['Permutation Intensity ',num2str(i_stim),' [.']);    
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        PermutationData                     = StimData(PERM(:,rep));          
        [Pval,StatVal,Coeffs]               = get_StatisticalValues(PermutationData,StimDesign);          
        Perm(i_stim).Pval(:,rep)            = Pval;        
        Perm(i_stim).Sval(:,rep)            = StatVal;
        Perm(i_stim).Coeffs(:,rep)          = Coeffs;
    end    
    utils.progressBar('1');   
    
    testCriterion                   = abs(squeeze(Test(i_stim).Coeffs));
    permCriterion                   = abs(squeeze(Perm(i_stim).Coeffs)); 
    P                               = crit2pval(testCriterion,permCriterion);
    Perm(i_stim).ContrastPval        = P;    
end
TestResults.PermutationPval = cat(1,Perm.ContrastPval);
end


if BOT_flag,
   disp('Starting permutation analysis!')
% BOOTSTRAP ANALYSIS 
Boot    = struct();
for i_stim=1:NUM_STIM,            
    StimData        = DATA(STIM==i_stim,:);
    StimDesign      = DESIGN_MATRIX(STIM==i_stim,:);
    BotDesign       = DESIGN(STIM==i_stim,:);
    BOT            	= get_BootMatrix(BotDesign,NUM_REP);
    
    utils.progressBar(['Bootstrap Intensity ',num2str(i_stim),' [.']);    
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        BootData                        = StimData(BOT(:,rep));        
        [Pval,StatVal,Coeffs,MargMeans] = get_StatisticalValues(BootData,StimDesign);          
        Boot(i_stim).Pval(:,rep)        = Pval;        
        Boot(i_stim).Sval(:,rep)        = StatVal;
        Boot(i_stim).Coeffs(:,rep)      = Coeffs;
        Boot(i_stim).MargMeans(:,rep)   = MargMeans;
    end
    utils.progressBar('1');           
    Boot(i_stim).Power       = nanmean(Boot(i_stim).Pval<ALPHA_ERROR,2);
    Boot(i_stim).CoeffsCI    = cat(2,quantile(Boot(i_stim).MargMeans,ALPHA_ERROR./2,2),quantile(Boot(i_stim).MargMeans,(1-ALPHA_ERROR./2),2));      
end
TestResults.Power       = single(cat(2,Boot.Power))';
TestResults.CoeffsCI    = permute(single(cat(3,Boot.CoeffsCI)),[3,1,2]);
end


%% ------------------------------------------------------------------------
% PERMUTATION BASED CLUSTER ANALYSIS 
% -------------------------------------------------------------------------
if CLUS_flag,
    disp('Starting Cluster analysis!')
% Reorder  into stim,factor,repetition
Perm_Pval = permute(cat(3,Perm.Pval),[3,1,2]); 
Perm_Sval = permute(cat(3,Perm.Sval),[3,1,2]);
for k=1:length(ClusterResults)
    TestClusterNum              = length(ClusterResults(k).Sval);    
    AbsTestClustVal             = abs(ClusterResults(k).Sval);
    
    if TestClusterNum==0,    
        ClusterResults(k).PermPval = [];
        disp(['Cluster Analysis Factor ',num2str(k),' is Empty!']);    
        continue; 
    end            
    
    if TestClusterNum>0  
                
        SigCounts                   = true(NUM_REP,TestClusterNum);
            
        utils.progressBar(['Cluster Analysis Factor ',num2str(k),' [.']);    
        for rep=1:NUM_REP,    
            utils.progressBar(rep);
            
            Hgrid               = Perm_Pval(:,k,rep)<ALPHA_ERROR;
            Sgrid               = Perm_Sval(:,k,rep);
            PermClusVal         = utils.stats2cluster(Hgrid,Sgrid);
            AbsPermClusVal      = sort(abs(PermClusVal));                        
            IsSignificantbyChance = true(1,TestClusterNum);
            
            if isempty(PermClusVal),
                IsSignificantbyChance   = false(1,TestClusterNum);
                SigCounts(rep,:)        = IsSignificantbyChance;  
                continue; 
            end

            for PickTestCluster=1:TestClusterNum
                IsSignificantbyChance(PickTestCluster) = AbsPermClusVal(1)>AbsTestClustVal(PickTestCluster);
            end
            SigCounts(rep,:)    = IsSignificantbyChance;              
        end   
        utils.progressBar('1');   
            
        tmpPval         = sum(SigCounts)./(NUM_REP);
        ClusterPermPval = max(tmpPval,1/NUM_REP);   
    end    
    ClusterResults(k).PermPval =  ClusterPermPval;  
end

end

end
%% LOCAL FUNCTIONS
% ANALYSIS FUNCTION 
function [Pval,StatVal,Coeffs,MargMeans] = get_StatisticalValues(DATA,DESIGN_MATRIX)   
    % for k=1:8, bin2dec(nms{k}([3,8,13]+1)), end -> [0 4 2 6 1 5 3 7]
    % designmatrix is [6 7 4 5 2 3 0 1]
    % f = []; for k=1:8, f = [f,find(ismember(a,b(k)))], end -> [4 8 2 6 3 7 1 5]; 
    RestoreOrder = [4 8 2 6 3 7 1 5];
    % DEF_MODEL will returning Interaction 
    % Interaction explained as Crosstable:
    %      90  45               90 45
    % M1    +  -     oder  M1   -   +
    % NPMA  -  +          NPMA  +   -
    % We show the sign of BI 90° and MO 45°   
    DEF_MODEL               = [1 0 0 0;0 1 0 0;0 0 1 0;1 0 1 0;0 1 1 0;0 0 0 1];        
    % Picking Values will be for 'Waveform', 'Orientation', 'Target','Interaction Waveform x Target',,'Interaction Orientation x Target')
    PickVal                 = 2:6;    
    % Picking Coeffs will be for ('BI=1', 'LM=1', 'M1 =1','BI=1 * M1=1' and 'LM=1 * M1=1')  % stats.coeffnames(PickCoeffs)
    PickCoeffs              = [3,5,7,11,15];    
    statistical_test        = @(x,y)anovan(x,y,'random',4,'display','off','varnames',{'BI','LM','M1','SUB'},'model',DEF_MODEL);        
    [~,tab,stats]           = statistical_test(DATA,DESIGN_MATRIX);
    
    Pval                    = single([tab{PickVal,7}]);
    StatVal                 = single([tab{PickVal,6}]);   
    Coeffs                  = single(stats.coeffs(PickCoeffs));

    warning('off','stats:multcompare:IgnoringRandomEffects')
    [~,m,~,~]               = multcompare(stats,'dim',[1,2,3],'display','off');    
    warning('on','stats:multcompare:IgnoringRandomEffects')
    MargMeans                = m(RestoreOrder,1); % in order of cat(1,setup.IO.BI,setup.IO.LM,setup.IO.M1)', i.e. DesignMatrix
end

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

% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Draws randomly from subjects with replacement for all factor levels present in the true data 
% @returns BOT, a factor-balanced bootstrap matrix
function BOT = get_BootMatrix(DESIGN,NUM_REP)
    BOT         = int32([]);
    dec_design  = bin2dec(num2str(DESIGN))+1;
    u_dsg       = unique(dec_design);
    for rep=1:NUM_REP, 
        bot_set    = NaN(size(dec_design));
        for i_dsg=1:length(u_dsg)
            dsg_GetPut  = dec_design==u_dsg(i_dsg);
            bot_sample  = datasample(find(dsg_GetPut),sum(dsg_GetPut));
            bot_set(dsg_GetPut) = bot_sample;
        end
        % output test: all(((dec_design(bot_set))-dec_design) == 0)
        BOT(:,rep) = bot_set;
    end
end


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
