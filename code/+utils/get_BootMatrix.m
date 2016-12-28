% CONSTRUCTION OF THE BOOTSTRAPPING VECTOR
% Draws randomly from subjects with replacement for all factor levels present in the true data 
% @returns BOT, a factor-balanced bootstrap matrix
function BOT = get_BootMatrix(DESIGN,NUM_REP)
    BOT         = int32([]);
    dec_design  = bin2dec(num2str(DESIGN))+1;
    u_dsg       = unique(dec_design);
    for rep=1:NUM_REP
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