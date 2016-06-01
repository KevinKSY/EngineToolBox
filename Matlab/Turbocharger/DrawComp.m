function DrawComp(SAE_map_comp,m_dot_user,pr_user,Level,axe)
    if isempty(axe)
        figure
        axe = axes;
    end;
    RPM = [0;SAE_map_comp.RPM];
    PR = [1;SAE_map_comp.PR];
    eff = [0;SAE_map_comp.eff];
    m_dot = [0;SAE_map_comp.m_dot];
    m_dot_surge = [0;SAE_map_comp.m_dot_surge];
    m_dot_chock = [0;SAE_map_comp.m_dot_chock];
    PR_surge = [1;SAE_map_comp.PR_surge];
    PR_chock = [1;SAE_map_comp.PR_chock];
    plot(axe,m_dot_surge,PR_surge,'r');
    hold on
    plot(axe,m_dot_chock,PR_chock,'r');
    hold on
    RPM_rep = unique(RPM);
    m = length(RPM_rep);

    for i=1:m
        idx = find(RPM == RPM_rep(i));
        plot(axe,m_dot(idx),PR(idx),'-*');
        hold on
    end;

    m_dot_min = 0;
    m_dot_max = max(m_dot);
    m_dot_fit = linspace(m_dot_min,m_dot_max,100);
    PR_min = 1.0;
    PR_max = max(PR);
    PR_fit = linspace(PR_min,PR_max,100);
    [eff_grid] = gridfit(m_dot,PR,eff,m_dot_fit,PR_fit);
    if isempty(Level)
        Level = [0.6:0.05:0.7 0.72:0.02:0.84 0.85];
    end;
    [~,h] = contour(axe,m_dot_fit,PR_fit,eff_grid,Level);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')/2) 
    hold on
    if ~isempty(m_dot_user)
        dmCorr = m_dot_user.*(1e5/pu).*sqrt(Tu/298);
        plot(axe,dmCorr,pr_user,'*');
    end;
    hold off
end