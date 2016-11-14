function [bot_CI,P] = do_botandperm_IOC(subAverage,DATA,PERM,BOT,do_test,num_rep)

% PERFORMING REAL ANALYSIS
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
true_VAL = [];
for si=1:7,
    t_sel               = (subAverage.Design(:,end)==si);
    t_design            = subAverage.Design(t_sel,1:4);      
    t_idx               = (subAverage.Design(:,end)==si);
    t_values            = DATA(t_idx,:);
    [p,tab,stats]       = do_test(t_values,t_design);
    true_VAL(si,:)      = [tab{2:5,6}];
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

% PERFORMING BOOTSTRAP ANALYSIS
% Output: bot_VAL(repetition,intensity,level,factor)
% values are the bootstrapped average for each level (true,false) of each factor
bot_VAL = [];
for si=1:7,
    for rep=1:num_rep,
        t_idx               = BOT{si}(:,rep);
        t_values            = DATA(t_idx);   
        t_design            = subAverage.Design(t_idx,1:4);
        a                   = [nanmean(t_values(t_design(:,1)==1)),nanmean(t_values(t_design(:,2)==1)),nanmean(t_values(t_design(:,3)==1))];
        b                   = [nanmean(t_values(t_design(:,1)==0)),nanmean(t_values(t_design(:,2)==0)),nanmean(t_values(t_design(:,3)==0))];
        bot_VAL(rep,si,:,:) = cat(1,a,b);
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

% ESTIMATION OF DESCRIPTIVE PARAMETERS
% Output: A matrix containing for each intensity the bootstrapped upper and lower  95% CI for each
% factor and level (CiLow,CiUp,StimulationIntensity,Level,Factor)
bot_CI = cat(1,quantile(bot_VAL,0.025,1),quantile(bot_VAL,0.975,1));