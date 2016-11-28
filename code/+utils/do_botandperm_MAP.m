function [true_M,true_P,perm_P,bot_CI,bot_P] = do_botandperm_MAP(DATA,SUBID,DESIGN,NUM_REP)

clc
if nargout<3,
    flag_nonparametric = false;
else
    flag_nonparametric = true;
end
data_dim        = size(DATA);

% DEFINITION OF ANALYSIS FUNCTION
NUM_REP         = 1000;
DEF_MODEL       = [1 0 0;0 1 0;1 1 0;0 0 1];
DESIGN_MATRIX   = cat(2,DESIGN,SUBID);
PICK_K          = 2:5;
do_test         = @(x,y)anovan(x,y,'random',3,'display','off','varnames',{'BI','LM','SUB'},'model',DEF_MODEL);
%% ------------------------------------------------------------------------
% PERFORMING TRUE ANALYSIS
%
% Output: true_VAL(intensity,factor)
% analysis function based on earlier anonymous definition
% Output: true_M(intensity,level,factor) based on marginal average

true_STATVAL    = single([]);
true_P          = single([]);
true_M          = single([]);
for i_pos=1:size(DATA,1),     
    t_data                  = (DATA(i_pos,:))';
    [~,tab,stats]           = do_test(t_data,DESIGN_MATRIX);
    true_STATVAL(i_pos,:)   = [tab{PICK_K,6}];
    true_P(i_pos,:)         = [tab{PICK_K,7}];    
    [~,m,~,~]               = multcompare(stats,'dim',[1,2],'display','off');
    true_M                  = cat(2,true_M,m(:,1));
end
true_M = true_M';

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
% Output: A cell array containing for each stimulation intensity a factor-balanced bootstrap matrix
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
% Output: perm_VAL(repetition,intensity,factor)
% analysis function based on earlier anonymous definition
perm_P          = ones(size(DATA,1),length(PICK_K));
for i_pos=1:size(DATA,1),
    perm_STATVAL = [];
    utils.progressBar('0');
    for rep=1:NUM_REP,         
        utils.progressBar(rep);
        t_data                  = (DATA(i_pos,PERM(:,rep)));            
        [~,tab,~]               = do_test(t_data,DESIGN_MATRIX);                
        perm_STATVAL            = cat(1,perm_STATVAL,[tab{PICK_K,6}]);
    end
    utils.progressBar('1');
    tmpP            = mean(repmat(true_STATVAL(i_pos,:),NUM_REP,1)>=perm_STATVAL);
    perm_P(i_pos,:) = tmpP;
end
%% ------------------------------------------------------------------------
% BOOTSTRAP ANALYSIS 
%
% PERFORMING
% Output: bot_VAL(repetition,intensity,level,factor) values are the 
% bootstrapped marginal average for each level (true,false) of each factor

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
    bot_CI(i_pos,:,:)   = cat(2,quantile(tmp_M,0.025,2),quantile(tmp_M,0.975,2));
  
end
% 
% ALPHA_VAL       = 0.05;
% D               = (pdist2(utils.get_DesignGrid,utils.get_DesignGrid));
% Adjacency       = D<6;        
% for k=1:size(bot_P,2)
%     for rep=1:num_rep,      
%         SigVal          = (bot_P(:,k,rep))<=ALPHA_VAL;
%         SigMat          = (double(SigVal)*double(SigVal)');
%         
%         
%         
%         [r,c]           = find(SigMat),
%         Z               = zeros(size(Adjacency)); 
%         Z(r(r~=c),c(c~=r)) = true;
%         
%     end
% end
