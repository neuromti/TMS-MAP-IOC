function ANT = get_ANT(setup,sub)
    ANT = [setup.SUB.pos{sub.subID}{sub.condition+1,5:7}];
    if any(isnan(ANT)) || isempty(ANT)
        ANT = setup.ANT(sub.subID,:);
    end
end