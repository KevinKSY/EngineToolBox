function DrawTurb(SAE_map_turb,m_dot_user,pr_user,eff_user,axe)
    if isempty(axe)
        figure
        axe = axes;
    end;
    RPM = [SAE_map_turb.RPM];
    PR = [SAE_map_turb.PR];
    eff = [SAE_map_turb.eff];
    m_dot = [SAE_map_turb.m_dot];
    RPM_rep = unique(RPM);
    m = length(RPM_rep);
    for i=1:m
        idx = find(RPM == RPM_rep(i));
        h = plotyy(PR(idx),m_dot(idx),PR(idx),eff(idx));
        set(h(2),'yLim',[0.7 0.9]);
        hold on
    end;
       
    %{
    m_dot_min = 0;
    m_dot_max = max(m_dot);
    m_dot_fit = linspace(m_dot_min,m_dot_max,100);
    PR_min = 1.0;
    PR_max = max(PR);
    PR_fit = linspace(PR_min,PR_max,100);
    [eff_grid] = gridfit(PR,m_dot,eff,PR_fit,m_dot_fit);
    if isempty(Level)
        Level = [0.6:0.05:0.7 0.72:0.02:0.84 0.85];
    end;
    [~,h] = contour(axe,PR_fit,m_dot_fit,eff_grid,Level);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')/2) 
    hold on
    %}
    if ~isempty(m_dot_user)
        [h,hLine1,hLine2] = plotyy(pr_user,m_dot_user,pr_user,eff_user);
        hLine1.LineStyle = 'none';
        hLine1.Marker = '*';
        hLine2.LineStyle = 'none';
        hLine2.Marker = 'o';
        set(h(2),'yLim',[0.7 0.9]);
    end;
    hold off
    set(h(1),'yLim',[floor(0.5*min(m_dot)) ceil(1.5*max(m_dot))]);
end