function [TestResults,ClusterResults] = mappingdata2statistics(DATA,SUBID,DESIGN,PERM_flag,BOOT_flag)
%% -------------------------------------------------------------------------
% DEFINITION OF PARAMETERS
% -------------------------------------------------------------------------
ALPHA_ERROR         = utils.get_ALPHAERROR();
DESIGN_MATRIX       = cat(2,DESIGN,SUBID);
NUM_REP             = 1000;
disp(['Alpha Error: ',num2str(ALPHA_ERROR)])
disp(['Number of Repetitions ',num2str(NUM_REP)])

% Flags for structured code testing 
if nargin <5,
    PERM_flag   = true;
    BOOT_flag   = true;
end
if nargout >1 && PERM_flag
    CLUS_flag   = true;
else
    CLUS_flag   = false;
end

%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
% -------------------------------------------------------------------------

for i_pos=1:size(DATA,1),     
    PositionData                    = (DATA(i_pos,:))';
    [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(PositionData,DESIGN_MATRIX);       
    Test(i_pos).Pval                = Pval;        
    Test(i_pos).Sval                = StatVal;
    Test(i_pos).Coeffs              = Coeffs';
    Test(i_pos).MargMeans           = MargMeans;
end

Test_Pval   = cat(1,Test.Pval);
Test_Sval   = cat(1,Test.Sval);
Test_Coeffs = cat(1,Test.Coeffs);
Test_MargMeans = cat(2,Test.MargMeans);

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

if PERM_flag
% PERMUTATION ANALYSIS 
PERM    = utils.get_PermMatrix(SUBID,NUM_REP);
Perm = struct();
for i_pos=1:size(DATA,1),        
   
  %  utils.progressBar(['Permutation GridPoint ',num2str(i_pos),' [.']);    
    for rep=1:NUM_REP,         
  %     utils.progressBar(rep);
        PermutationData                     = (DATA(i_pos,PERM(:,rep)))';          
        [Pval,StatVal,MargMeans,Coeffs]     = get_StatisticalValues(PermutationData,DESIGN_MATRIX);          
        Perm(i_pos).Pval(:,rep)             = Pval;        
        Perm(i_pos).Sval(:,rep)             = StatVal;
        Perm(i_pos).Coeffs(:,rep)           = Coeffs;
        Perm(i_pos).MargMeans(:,rep)        = MargMeans;
    end    
  %  utils.progressBar('1');   
    
    testCriterion                   = abs(squeeze(Test(i_pos).Coeffs));
    permCriterion                   = abs(squeeze(Perm(i_pos).Coeffs)); 
    P                               = utils.crit2pval(testCriterion,permCriterion);
    Perm(i_pos).ContrastPval        = P;    
end
TestResults.PermutationPval = cat(1,Perm.ContrastPval);
end


if BOOT_flag,
% BOOTSTRAP ANALYSIS 
BOT     = utils.get_BootMatrix(DESIGN,NUM_REP);
Boot    = struct();
for i_pos=1:size(DATA,1),                
 %   utils.progressBar(['Bootstrap GridPoint ',num2str(i_pos),' [.']);    
    for rep=1:NUM_REP,         
 %       utils.progressBar(rep);
        BootData                        = (DATA(i_pos,BOT(:,rep)))';
        BootDesign                      = cat(2,DESIGN,SUBID(BOT(:,rep)));
        [Pval,StatVal,MargMeans,Coeffs] = get_StatisticalValues(BootData,BootDesign);          
        Boot(i_pos).Pval(:,rep)          = Pval;        
        Boot(i_pos).Sval(:,rep)          = StatVal;
        Boot(i_pos).Coeffs(:,rep)        = Coeffs;
        Boot(i_pos).MargMeans(:,rep)     = MargMeans;
    end
 %   utils.progressBar('1');           
    Boot(i_pos).Power       = nanmean(Boot(i_pos).Pval<ALPHA_ERROR,2);
    Boot(i_pos).CoeffsCI    = cat(2,quantile(Boot(i_pos).Coeffs,ALPHA_ERROR./2,2),quantile(Boot(i_pos).Coeffs,(1-ALPHA_ERROR./2),2));      
end
TestResults.Power       = single(cat(2,Boot.Power));
TestResults.CI          = single(cat(3,Boot.CoeffsCI));
end
%% ------------------------------------------------------------------------
% PERMUTATION BASED CLUSTER ANALYSIS 
% -------------------------------------------------------------------------
if CLUS_flag == true
    
Perm_Pval = permute(cat(3,Perm.Pval),[3,1,2]);
Perm_Sval = permute(cat(3,Perm.Sval),[3,1,2]);
for k = 1 : length(ClusterResults)
    TestClusterNum              = length(ClusterResults(k).Sval);    
    AbsTestClustVal             = abs(ClusterResults(k).Sval);
    if TestClusterNum == 0    
        ClusterResults(k).PermPval = [];
        disp(['Cluster Analysis Factor ',num2str(k),' is Empty!']);    
        continue; 
    end            
    
    if TestClusterNum>0  
                
        SigCounts                   = true(NUM_REP,TestClusterNum);
            
  %      utils.progressBar(['Cluster Analysis Factor ',num2str(k),' [.']);    
        for rep = 1 : NUM_REP  
  %         utils.progressBar(rep);
            
            Hgrid               = utils.vec2mesh(Perm_Pval(:,k,rep)<ALPHA_ERROR);
            Sgrid               = utils.vec2mesh(Perm_Sval(:,k,rep));
            PermClusVal         = utils.stats2cluster(Hgrid,Sgrid);
            AbsPermClusVal   = sort(abs(PermClusVal));                        
            IsSignificantbyChance = true(1,TestClusterNum);
            
            if isempty(PermClusVal)
                IsSignificantbyChance   = false(1,TestClusterNum);
                SigCounts(rep,:)        = IsSignificantbyChance;  
                continue; 
            end

            % Compare against highest value
            PickPermCluster = 1;
            for PickTestCluster = 1 : TestClusterNum
                IsSignificantbyChance(PickTestCluster) = AbsPermClusVal(PickPermCluster) > AbsTestClustVal(PickTestCluster);
                if AbsPermClusVal(PickPermCluster) > AbsTestClustVal(PickTestCluster)
                    PickPermCluster = min(PickPermCluster+1,length(AbsPermClusVal));
                end
            end
            SigCounts(rep,:)    = IsSignificantbyChance;              
        end   
   %     utils.progressBar('1');   
            
        tmpPval         = sum(SigCounts)./(NUM_REP);
        ClusterPermPval = max(tmpPval,1/NUM_REP);   
    end    
    ClusterResults(k).PermPval =  ClusterPermPval;  
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