function [bot_CI,P,true_M] = do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep)
%supress warning that random effects are ignored when estimating the means
warning('off','stats:multcompare:IgnoringRandomEffects');
% PERFORMING REAL ANALYSIS
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
% Output: true_M(intensity,level,factor) based on marginal average
true_VAL = [];
true_M  = [];
for si=1:7,
    t_sel               = (subAverage.Design(:,end)==si);
    t_design            = subAverage.Design(t_sel,1:4);      
    t_idx               = (subAverage.Design(:,end)==si);
    t_values            = DATA(t_idx,:);
    [~,tab,stats]       = do_test(t_values,t_design);
    [~,m1]              = multcompare(stats,'dim',1,'display','off');
    [~,m2]              = multcompare(stats,'dim',2,'display','off');
    [~,m3]              = multcompare(stats,'dim',3,'display','off');
    b                   = [m1(1,1),m2(1,1),m3(1,1)]; %factor level: false 
    a                   = [m1(2,1),m2(2,1),m3(2,1)]; %factor level: true
    true_M(si,:,:)      = cat(1,a,b);
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
bot_VAL = [];
for si=1:7,
    for rep=1:num_rep,
        t_idx               = BOT{si}(:,rep);
        t_values            = DATA(t_idx);   
        t_design            = subAverage.Design(t_idx,1:4);
        
        [~,~,stats]         = do_test(t_values,t_design);
        [~,m1]              = multcompare(stats,'dim',1,'display','off');
        [~,m2]              = multcompare(stats,'dim',2,'display','off');
        [~,m3]              = multcompare(stats,'dim',3,'display','off');
        b                   = [m1(1,1),m2(1,1),m3(1,1)]; %factor level: false 
        a                   = [m1(2,1),m2(2,1),m3(2,1)]; %factor level: true
        bot_VAL(rep,si,:,:) = cat(1,a,b);
    end   
end

% ESTIMATION OF DESCRIPTIVE PARAMETERS
% Output: A matrix containing for each intensity the bootstrapped upper 
% and lower  95% CI for each factor and level 
% (CiLow,CiUp,StimulationIntensity,Level,Factor)
bot_CI = cat(1,quantile(bot_VAL,0.025,1),quantile(bot_VAL,0.975,1));