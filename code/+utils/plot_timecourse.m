function h = plot_timecourse(MargMeans,fig_DESIGN,fig_LABEL,Pval)

    % Assert Parameters are valid
    if ~islogical(fig_DESIGN), warning('pltTC:ForceLogical','Design Matrix not logical. Forcing Type Change'), fig_DESIGN = logical(fig_DESIGN); end
    assert(size(fig_DESIGN,1) == 2, 'Design Matrix wrong dimensions');  
    assert(iscell(fig_LABEL),'Label not a Cell Array'); 
    assert(size(fig_LABEL,2) == 3,'Label Cell Array wrong dimensions');

    
    
    % Aggregate
    
    ave1    = squeeze(mean(MargMeans(:,fig_DESIGN(1,:),:),2));  
    ave2    = squeeze(mean(MargMeans(:,fig_DESIGN(2,:),:),2));

    PosVal = cat(2,ave1,ave2);
   
    % Plot
    
    clear h
    h(1) = plot_tc(ave1,fig_LABEL{1});
    h(2) = plot_tc(ave2,fig_LABEL{2});
    h(3) = plot_sig(ave1,ave2,Pval,fig_LABEL);
    
end

function h = plot_tc(valu,tit)
   
    CLIM        = -125:5:125;    
  
    h = plot_whatever(valu,CLIM,tit,14);
    
end

function h = plot_sig(valu1,valu2,Pval,tit)
    
    valu            = (valu1-valu2);
    sig             = -log10(flipud(Pval)).*sign(valu);
    CLIM           = [-3:0.01:-1.3,0,1.3:0.01:3];
    tit             = sprintf('%s vs %s ',tit{1:2});
    h = plot_whatever(sig,CLIM,tit,14);
        
end

function h = plot_whatever(PlotValu,CLIM,tit,FontSize)
    PlotValu       = cat(2,PlotValu,linspace(CLIM(1),CLIM(end),7)');
    
    h = figure;
    hold on
    set(gcf,'Position',[100,100,1000,450],'paperpositionmode','auto')
    [~,ch] = contour(PlotValu,repmat(CLIM,1,2),'LevelListMode','manual');
    set(ch,'fill','on')
        
    set(gca,'XLIM',[0 275],'XTICK',0:25:275,'XTICKLABEL',5:5:60,'fontsize',FontSize)
    set(gca,'YLIM',[1 7],'YTICK',1:7,'YTICKLABEL',90:10:150,'fontsize',FontSize)
    ylabel('Stimulation Intensity in % RMT','fontsize',FontSize)
    xlabel('Latency in ms','fontsize',FontSize)
    title(tit,'fontsize',FontSize)   
end