%> @function
%> Plots two input-output curves in two different colors,
%> Plots patch for confidence intervals
%> @params 

function plot_ioc(MargMeans,CI,sel_DESIGN,Labels,Pval)

    % Assert Parameters are valid
    if ~islogical(sel_DESIGN), warning('pltIOC:ForceLogic','Design Matrix not logical. Forcing Type Change'), sel_DESIGN = logical(sel_DESIGN); end
    assert(size(sel_DESIGN,1) == 2, 'Design Matrix wrong dimensions');  
    assert(iscell(Labels),'Label not a Cell Array'); 
    assert(size(Labels,2) == 3,'Label Cell Array wrong dimensions');

    
    % Plot
    figure
    set(gcf,'Position',[100,100,800,450],'paperpositionmode','auto')
    hold on
    
    PosVal  = cat(2,mean(MargMeans(:,sel_DESIGN(1,:)),2),mean(MargMeans(:,sel_DESIGN(2,:)),2));
    ave     = mean(MargMeans(:,sel_DESIGN(1,:)),2);
    ciup    = mean(CI(:,sel_DESIGN(1,:),2),2);
    cilo    = mean(CI(:,sel_DESIGN(1,:),1),2);
    h_a1    = plot_single_ioc(ave,ciup,cilo,'r');
    
    ave     = mean(MargMeans(:,sel_DESIGN(2,:)),2);
    ciup    = mean(CI(:,sel_DESIGN(2,:),2),2);
    cilo    = mean(CI(:,sel_DESIGN(2,:),1),2);
    h_a2    = plot_single_ioc(ave,ciup,cilo,'b');

    % Annotate / Beautify    
    title([[Labels{1},' vs. ',Labels{2}],' -> ',Labels{3}])
    legend([h_a1,h_a2],Labels(1:2),'location','northwest')    
    
    if any(regexpi(Labels{3},'Amplitude'))
        set(gca,'YLIM',[0 700],'YTICK',0:100:700,'YTICKLABEL',0:100:700)
        ylabel('Amplitude in µVpp')
    elseif any(regexpi(Labels{3},'MEP'))
        ylabel('MEP in %')            
        set(gca,'YLIM',[0 1],'YTICK',0:0.1:1,'YTICKLABEL',0:.1:1)
    elseif any(regexpi(Labels{3},'Latency'))
        ylabel('Latency in ms')
        set(gca,'YLIM',[24.5 28.5],'YTICK',25:1:28,'YTICKLABEL',25:1:28)         
    else
        warning('pltIOC:UnitUnrecognized','No recognized unit of measure')            
    end
    
    set(gca,'XLIM',[0.5 7.5],'XTICK',1:7,'XTICKLABEL',90:10:150,'fontsize',14)
    xlabel('Stimulation Intensity in % RMT')
    plot_sigstars(Pval,PosVal);
    
end
%%-------------------------------------------------------------------------
% plots single input-output curve
function [h_a,h_p] = plot_single_ioc(ave,ciup,cilo,fig_color)
    
    h_p     = patch([1:7,7:-1:1],[cilo',flipud(ciup)'],ones(1,14));
    h_a     = plot(ave);
    set(h_a,'Linewidth',2,'color',fig_color,'marker','s','markeredgecolor','k','markersize',3)
    set(h_p,'Linewidth',1,'facecolor',fig_color,'facealpha',0.15,'edgealpha',0)
end

function plot_sigstars(Pval,PosVal)
    SigIdx = Pval < utils.get_ALPHAERROR(length(Pval));
    if any(SigIdx)
        PlotIdx  = find(SigIdx);
        if size(PosVal,1) == 1
            plot(PlotIdx,PosVal,'k*','markerfacecolor','k','markersize',10);
        elseif size(PosVal,1) == size(Pval,1)   
            if size(PosVal,2) == 2
                plot(PlotIdx,PosVal(PlotIdx,1),'linestyle','none','marker','o','markerfacecolor','r','markeredgecolor','r','markersize',10);
                plot(PlotIdx,PosVal(PlotIdx,1),'ks','markerfacecolor','k','markersize',4);
                plot(PlotIdx,PosVal(PlotIdx,2),'linestyle','none','marker','o','markerfacecolor','b','markeredgecolor','b','markersize',10);
                plot(PlotIdx,PosVal(PlotIdx,2),'ks','markerfacecolor','k','markersize',4);                
            elseif size(PosVal,2) == 1
                plot(PlotIdx,PosVal(PlotIdx,1),'k*','markerfacecolor','k','markersize',10);
                
            else
                warning('pltIOC:WrongDim','Wrong Dimensions')
            end
        else
            warning('pltIOC:WrongDim','Wrong Dimensions')
        end       
    end        
end
