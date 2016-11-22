function [bot_CI,bot_D,P,true_M] = do_botandperm_MAP(DATA,SUB,DESIGN,num_rep,do_test)

data_dim    = size(DATA);
DATA        = reshape(DATA,prod(data_dim(1:3)),[]);

% CONSTRUCTION OF THE PERMUTATION ARRAY
% Output: A cell array containing for each stimulation intensity a matrix for later design matrix permutation 
PERM        = [];
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
BOT         = {};
dec_design  = bin2dec(num2str(DESIGN))+1;
u_dsg       = unique(dec_design);
for rep=1:num_rep, 
    bot_set    = [];
    for i_dsg=1:length(u_dsg)
        dsg_GetPut  = dec_design==u_dsg(i_dsg);
        bot_set     = cat(1,bot_set,datasample(find(dsg_GetPut),sum(dsg_GetPut)));
    end
    BOT{rep} = bot_set;
end

% PERFORMING TRUE ANALYSIS
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
% Output: true_M(intensity,level,factor) based on marginal average

TODO :MAKE FASTER!

%--------------------------------------------------------------------------
% PERFORMING PERMUTATION ANALYSIS
% Output: perm_VAL(repetition,intensity,factor)
% analysis function based on earlier anonymous definition
t_design = cat(2,DESIGN);
perm_M   = [];
for i_pos=1:size(DATA,1),
    t_values            = DATA(i_pos,:);  
    t_values            = (t_values(PERM));
    do_grp              = @(x)[grpstats(x,t_design(:,1));grpstats(x,t_design(:,2))]';
    m                   = (do_grp(t_values));
    perm_M              = cat(3,perm_M,cat(2,m(:,1)-m(:,2),m(:,3)-m(:,4)));  
end











true_VAL    = [];
true_M      = [];
for si=1:7,
    t_sel               = (subAverage.Design(:,end)==si);
    t_design            = subAverage.Design(t_sel,1:4);      
    t_idx               = (subAverage.Design(:,end)==si);
    t_values            = DATA(t_idx,:);
    [~,tab,~]           = do_test(t_values,t_design);
    

    true_M(si,:,:)      = marginal_m;
    
    true_VAL(si,:)      = [tab{2:5,6}];
end

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% PERFORMING BOOTSTRAP ANALYSIS
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
