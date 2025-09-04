function p = plot_finer_coh_grid(D, coh_fine, Pfine, ndt_m)

    p = publish_plot(2,1);
    set(gcf,'Position',[263   134   463   563]);

    choice = D.choice;
    coh = D.coh;
    rt = D.rt;
    c = D.c;
    

    p.next();

    [tt,xx,ss] = curva_media(choice,coh,[],0);
    errorbar(tt,xx,ss,'LineStyle','none','Marker','.','Color','k','MarkerSize',20);
    hold all

    plot(coh_fine,Pfine.up.p,'k-');
    xlim([min(coh_fine),max(coh_fine)]);
    xlabel('Motion coherence');
    ylabel('P rightward choice')

    p.next();

    % only correct trials
    rt_model_c = Pfine.up.mean_t;
    rt_model_c(coh_fine<0) = Pfine.lo.mean_t(coh_fine<0);
    rt_model_c(coh_fine==0) = (Pfine.up.mean_t(coh_fine==0)+Pfine.lo.mean_t(coh_fine==0))/2;
    rt_model_c = rt_model_c + ndt_m;

    [tt,xx,ss] = curva_media(rt,coh,c==1,0);
    errorbar(tt,xx,ss,'LineStyle','none','Marker','.','Color','k','MarkerSize',20);
    hold all
    plot(coh_fine,rt_model_c,'k-');
    xlim([min(coh_fine),max(coh_fine)]);
    xlabel('Motion coherence');
    ylabel('RT (s)')

    %set(gcf,'Position',[270   793  1084   293])

    p.format();

    
end
