function HS = get_HS(setup,sub)

    HS = [setup.SUB.pos{sub.subID}{sub.condition+1,2:4}];
   if any(isnan(HS)) || isempty(HS)
       HS = setup.HS(sub.subID,:);
   end
end