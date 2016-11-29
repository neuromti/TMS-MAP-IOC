function [parametric,permutated,bootstrapped] = mappingdata2statistics(DATA,SUBID,DESIGN)

clc
if nargout<2,
    flag_nonparametric = false;
else
    flag_nonparametric = true;
end

% DEFINITION OF ANALYSIS FUNCTION
NUM_REP         = 1000;
DEF_MODEL       = [1 0 0;0 1 0;1 1 0;0 0 1];
DESIGN_MATRIX   = cat(2,DESIGN,SUBID);
PICK_K          = 2:4;
PICK_S          = [3,5,9]; %Supposed to be ('BI=1', 'LM=1', 'BI=1 * LM=1') Literals instead of String Indexing was used for Speed
%PICK_S          = {'BI=1', 'LM=1', 'BI=1 * LM=1'}
% Interaction Explained as Crosstable:
%     90  45              90 45
% BI   +  -     oder  BI  -   +
% MO   -  +           MO  +   -
% We show the sign of BI 90° and MO 45°

ALPHA_ERROR     = 0.05;
do_test         = @(x,y)anovan(x,y,'random',3,'display','off','varnames',{'BI','LM','SUB'},'model',DEF_MODEL);
%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
%
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
% Output: true_M(intensity,level,factor) based on marginal average
warning('off')
true_STATVAL    = single([]);
true_P          = single([]);
true_S          = single([]);
true_M          = single([]);
for i_pos=1:size(DATA,1),     
    t_data                  = (DATA(i_pos,:))';
    [~,tab,stats]           = do_test(t_data,DESIGN_MATRIX);
    true_STATVAL(i_pos,:)   = [tab{PICK_K,6}];
    true_P(i_pos,:)         = [tab{PICK_K,7}];        
    true_S(i_pos,:)         = stats.coeffs(PICK_S); %stats.coeffs(ismember(stats.coeffnames,PICK_S))
    [~,m,~,~]               = multcompare(stats,'dim',[1,2],'display','off');
    true_M                  = cat(2,true_M,m(:,1));
end
true_M = true_M';

true_ClusVal = {};
true_ClusIdx = {};
for k=1:size(true_S,2),    
    Hgrid               = utils.vec2mesh((true_P(:,k)<ALPHA_ERROR));
    Sgrid               = utils.vec2mesh(true_S(:,k));
    
    posH                = Hgrid&(Sgrid>0);
    negH                = Hgrid&(Sgrid<0);
    
    [PosClusVal,posGridIdx] = utils.stats2cluster(posH,Sgrid);
    [NegClusVal,negGridIdx] = utils.stats2cluster(negH,-Sgrid);
    true_ClusVal{k}            = {PosClusVal,NegClusVal};
    true_ClusIdx{k}         = cat(3,posGridIdx,negGridIdx);
    
end
       
parametric = struct('MargMeans',true_M,'P',true_P,'Coeff',true_S);
parametric.ClusterVal = true_ClusVal;
parametric.ClusterIdx = true_ClusIdx;
warning('on')
%%
if ~flag_nonparametric, 
    warning('You chose not to run permutation and bootstrap analysis!')
    return;
else
    warning('Starting permutation and bootstrap analysis!')
end
%% ------------------------------------------------------------------------
% MIXING MATRICES
%
% CONSTRUCTION OF THE PERMUTATION ARRAY
% Permutes the factor levels for all subjects present in the true data 
% @returns PERM, a subject-balanced permutation matrix 
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

% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Draws randomly from subjects with replacement for all factor levels present in the true data 
% @returns BOT, a factor-balanced bootstrap matrix
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

clear u_dsg i_dsg rep i_sub dsg_GetPut sub_GetPut perm_set bot_set dec_design u_sub
%% ------------------------------------------------------------------------
% PERMUTATION ANALYSIS 
%
% PERFORMING PERMUTATION ANALYSIS
% analysis function based on earlier defined anonymous definition do_test
% @returns Pval, an alpha-error estimate for each grid point
perm_P          = single(ones(size(DATA,1),length(PICK_K)));
perm_ClusterP      = single(ones(size(DATA,1),length(PICK_K),NUM_REP));
perm_STATVAL    = single(ones(size(DATA,1),length(PICK_K),NUM_REP));
for i_pos=1:size(DATA,1),        
    utils.progressBar('0');
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        t_data                  = (DATA(i_pos,PERM(:,rep)));            
        [p,tab,~]               = do_test(t_data,DESIGN_MATRIX);                
        perm_STATVAL(i_pos,:,rep)   = [tab{PICK_K,6}];
        perm_ClusterP(i_pos,:,rep)  = [tab{PICK_K,7}];
    end    
    utils.progressBar('1');
    tmpP            = mean(repmat(true_STATVAL(i_pos,:),NUM_REP,1)>=squeeze(perm_STATVAL(i_pos,:,:))');
    perm_P(i_pos,:) = tmpP;
end
permutated = struct('Pval',perm_P);


% PERFORMING CLUSTER ANALYSIS
% analysis function based on earlier estimated test statistics
% @returns ClusterPval, an alpha-error estimate for each cluster
final_ClusterP = {};
for k=1:length(PICK_K)
    
    TrueClusVal     = [parametric.ClusterVal{k}{1},parametric.ClusterVal{k}{2}];

    if ~isempty(TrueClusVal)
        SigCounts       = false(NUM_REP,length(TrueClusVal));        
        for rep=1:NUM_REP,        
            Hgrid           = utils.vec2mesh((perm_ClusterP(:,k,rep)<ALPHA_ERROR));
            Sgrid           = utils.vec2mesh(perm_STATVAL(:,k,rep));            
            PosClusVal      = utils.stats2cluster(Hgrid&(Sgrid>0),Sgrid);
            NegClusVal      = utils.stats2cluster(Hgrid&(Sgrid<0),-Sgrid);
            PermClusVal     = sort([PosClusVal,-NegClusVal],'descend');
      
            c_c             = 1;
            IsSignificant 	= [];
            for i_c = 1:length(TrueClusVal),
                if c_c<=length(PermClusVal),
                    if TrueClusVal(i_c)>PermClusVal(c_c),
                        c_c             = c_c+1;
                        IsSignificant   = [IsSignificant,i_c];
                    end
                else
                    IsSignificant   = [IsSignificant,i_c];
                end
            end
            SigCounts(rep,IsSignificant) = true;
        end    
        final_ClusterP{k} = 1-((sum(SigCounts)-1)./NUM_REP);
    else
        final_ClusterP{k} = [];
    end
end
% 
permutated.ClusterPval = final_ClusterP;
%% ------------------------------------------------------------------------
% BOOTSTRAP ANALYSIS 
%
% PERFORMING
% analysis function based on earlier defined anonymous definition do_test
% @returns bot_P, a power estimate (1-beta)
% @returns bot_CI, an estimate of the 95% confidence interval
bot_CI = NaN(size(DATA,1),length(PICK_K),2);
bot_P  = NaN(size(DATA,1),length(PICK_K),NUM_REP);
for i_pos=1:size(DATA,1),
    tmp_M = [];
    tmp_P = [];        
    utils.progressBar('0');
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        t_data                  = (DATA(i_pos,BOT(:,rep)));            
        [p,~,stats]             = do_test(t_data,DESIGN_MATRIX);                
        [~,m,~,~]               = multcompare(stats,'dim',[1,2],'display','off');
        tmp_M                   = cat(2,tmp_M,m(:,1));
        tmp_P                   = cat(2,tmp_P,p);        
    end
    utils.progressBar('1');
    bot_P(i_pos,:,:)    = tmp_P;    
    bot_CI(i_pos,:,:)   = cat(2,quantile(tmp_M,ALPHA_ERROR./2,2),quantile(tmp_M,(1-ALPHA_ERROR./2),2));  
end

bootstrapped = struct('P',bot_P,'CI',bot_CI);