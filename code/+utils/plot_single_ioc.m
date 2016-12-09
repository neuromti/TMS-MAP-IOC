function h_a = plot_single_ioc(ave,fig_color)

    h_a     = plot(ave);
    set(h_a,'Linewidth',2,'color',fig_color,'marker','s','markeredgecolor','k','markersize',3)
    
end