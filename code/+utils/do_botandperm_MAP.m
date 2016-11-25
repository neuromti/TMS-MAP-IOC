function [bot_CI,bot_D,P,true_MM] = do_botandperm_MAP(DATA,SUB,DESIGN,lcl_folder,num_rep)

flag_nonparametric = false;

if ~exist(lcl_folder,'dir'), mkdir(lcl_folder); end

data_dim        = size(DATA);

% DEFINITION OF ANALYSIS FUNCTION
% output values are f-values from anova
% subject as random factor
% no interactions
% independently for each SI due to assumption of heteroscedasticity
%do_test         = @(x,y)anovan(x,y,'random',3,'display','off','varnames',{'BI','LM','SUB'});
DEF_MODEL       = [1 0 0;0 1 0;1 1 0;0 0 1];
do_test         = @(x,y)anovan(x,y,'random',3,'display','off','varnames',{'BI','LM','SUB'},'model',DEF_MODEL);
% Added do_margmean for speed in estimation of marginal means w/o interactions
do_margmean     = @(x,y)[grpstats(x,y(:,1));grpstats(x,y(:,2))]';
%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
%
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
% Output: true_M(intensity,level,factor) based on marginal average

if ~exist([lcl_folder,'map_true_stats.mat'],'file')
    t_design        = cat(2,DESIGN,SUB);
    true_STATVAL    = single([]);
    true_MM         = single([]);
    true_P          = single([]);
    delete(gcp('nocreate'));
    obj_pool    = parpool(4);
    parfor i_pos=1:size(DATA,1),     
        %tic
        t_data                  = (DATA(i_pos,:))';
        [~,tab,~]               = do_test(t_data,t_design);
        true_MM(i_pos,:)        = do_margmean(t_data,t_design);    
        true_STATVAL(i_pos,:)   = [tab{2:5,6}];
        true_P(i_pos,:)         = [tab{2:5,7}];
        %utils.toc_to_estimate(toc,data_dim(1));    
    end
    delete(obj_pool);
    save([lcl_folder,'map_true_stats.mat'],'true_MM','true_STATVAL','true_P');
else
    load([lcl_folder,'map_true_stats.mat'],'true_MM','true_STATVAL','true_P');
end


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
% Output: A cell array containing for each stimulation intensity a matrix for later design matrix permutation 
PERM        = int8([]);
u_sub       = unique(SUB);
for rep=1:num_rep, 
    perm_set    = [];
    for i_sub=1:length(u_sub)
        sub_GetPut  = find(SUB==u_sub(i_sub));        
        perm_set    = cat(1,perm_set,sub_GetPut(randperm(length(sub_GetPut)),:));
    end
    PERM(:,rep) = perm_set;
end

% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Output: A cell array containing for each stimulation intensity a factor-balanced bootstrap matrix
BOT         = int8([]);
dec_design  = bin2dec(num2str(DESIGN))+1;
u_dsg       = unique(dec_design);
for rep=1:num_rep, 
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
% Output: perm_VAL(repetition,intensity,factor)
% analysis function based on earlier anonymous definition

if ~exist([lcl_folder,'map_perm_stats.mat'],'file')
    t_design        = cat(2,DESIGN,SUB);
    modelOrder      = int8(size(DEF_MODEL,1));
    tab_sel         = 2:modelOrder+1;
    perm_P          = ones(size(DATA,1),modelOrder);
    parfor i_pos=1:size(DATA,1),
        perm_STATVAL = zeros(1,4);
        for rep=1:num_rep,                      
            t_data                  = (DATA(i_pos,PERM(:,rep)));            
            [~,tab,~]               = do_test(t_data,t_design);            
            tmp                     = (true_STATVAL(i_pos,:,:) >= [tab{tab_sel,6}]);      
            perm_STATVAL            = perm_STATVAL+tmp;
        end
        perm_P(i_pos,:) = (perm_STATVAL./(num_rep+1));        
    end
    delete(obj_pool);
    save([lcl_folder,'map_perm_stats.mat'],'perm_P')
else
    load([lcl_folder,'map_perm_stats.mat'],'perm_P')
end


% PERFORMING PERMUTATION ANALYSIS
% Output: perm_VAL(repetition,intensity,factor)
% analysis function based on earlier anonymous definition
perm_VAL = [];
for si=1:7,
    t_sel               = subAverage.Design(:,end)==si;
    t_design            = subAverage.Design(t_sel,1:4);        
    for rep=1:num_rep,
        t_idx               = PERM{si}(:,rep);
        t_values            = DATA(t_idx);                         
        [p,tab,stats]       = do_test(t_values,t_design);
        perm_VAL(rep,si,:)  = [tab{2:5,6}];
    end   
end

% ESTIMATION OF P-VALUE
% Output: A matrix containing for each stimulation intensity the p-value of
% a two-sided hypothesis test
P = [];
for si=1:7,
    matrix_P    = repmat(true_VAL(si,:),num_rep,1)>=squeeze(perm_VAL(:,si,:));
    t_p         = mean(matrix_P);
    t_p         = min(t_p,1-t_p).*2;
    P(si,:)     = t_p;
end

%% ------------------------------------------------------------------------
% BOOTSTRAP ANALYSIS 
%
% PERFORMING
% Output: bot_VAL(repetition,intensity,level,factor) values are the 
% bootstrapped marginal average for each level (true,false) of each factor
bot_M = [];
for si=1:7,
    for rep=1:num_rep,
        t_idx               = BOT{si}(:,rep);
        t_values            = DATA(t_idx);   
        t_design            = subAverage.Design(t_idx,1:4);
        
        m                   = (grpstats(t_values,t_design(:,1:3)));
        marginal_m          = cat(2,grpstats(m,setup.IO.BI),grpstats(m,setup.IO.LM),grpstats(m,setup.IO.M1));   
 
        bot_M(rep,si,:,:) = marginal_m;
    end   
end

% ESTIMATION OF DESCRIPTIVE PARAMETERS
% Output: A matrix containing for each intensity the bootstrapped upper 
% and lower  95% CI for each factor and level 
% (CiLow,CiUp,StimulationIntensity,Level,Factor)
bot_CI          = cat(1,quantile(bot_M,0.025,1),quantile(bot_M,0.975,1));
bot_DELTA       = squeeze(diff(bot_M,[],3));
bot_D           = squeeze(cat(1,quantile(bot_DELTA,0.025,1),quantile(bot_DELTA,0.975,1)));
