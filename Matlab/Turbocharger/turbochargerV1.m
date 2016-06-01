
classdef turbochargerV1
    properties (Hidden)
        ModelType
        effComp_grid_fit;
        m_dotComp_grid_fit;
        effTurb_grid_fit;
        m_dotTurb_grid_fit;
        Grid_no = 100;
        RPMComp_fit;    
        RPMTurb_fit;
        PRComp_fit;
        PRTurb_fit;
        m_dotComp_fit;
        m_dotTurb_fit;
        coeff_m_dot_turb;
        K_pTF_coeff = [0.00378955656268277,1.43487448145079;-0.00905287405591890,-0.110762167507735;0.000676269940804610,-0.0338912838653670;0.00557633083751460,0.0144999652821938;-0.000716126707482597,-0.0178991462433519;-2.48117033240293e-05,0.00697027168505133];
        Cp_pTF_coeff = [3.75627710318873,952.182081186873;-9.70682873742566,152.561091919067;1.29744990668128,40.1088939509955;6.35935024827548,31.1378807865909;-1.51207506730126,94.3052521819761;0.0299339569666709,-6.22325708679949];
    end
    properties (SetAccess = private)
        SAE_map_comp;
        SAE_map_turb;
        TC_map;
        Uc_opt;
    end
    methods
        function SetUc_opt(TC,Value)
            TC.Uc_opt = Value;
        end
        function Value = GetUc_opt(TC)
            Value = TC.Uc_opt;
        end
        function TC = turbochargerV1(SAE_map_comp,SAE_map_turb,compRef,turbRef)
            TC.SAE_map_comp = SAE_map_comp;
            TC.SAE_map_turb = SAE_map_turb;
           %% Compressor Calculation
            [~,TC.SAE_map_comp] = CompMapGenerator(TC.SAE_map_comp,compRef);
            m_dot = [0;TC.SAE_map_comp.m_dot];
            PR = [1;TC.SAE_map_comp.PR];
            RPM = [0;TC.SAE_map_comp.RPM];
            eff = [0;TC.SAE_map_comp.eff];
            m_dot_min = 0;
            m_dot_max = max(m_dot);
            TC.m_dotComp_fit = linspace(m_dot_min,m_dot_max,TC.Grid_no+1);
            PR_min = 1.01;
            PR_max = max(PR);
            TC.PRComp_fit = linspace(PR_min,PR_max,TC.Grid_no+1);
            RPM_min = 0;
            RPM_max = max(RPM);
            TC.RPMComp_fit = linspace(RPM_min,RPM_max,TC.Grid_no+1);
            [TC.m_dotComp_grid_fit] = gridfit(RPM,PR,m_dot,TC.RPMComp_fit,TC.PRComp_fit);
            [TC.effComp_grid_fit] = gridfit(m_dot,PR,eff,TC.m_dotComp_fit,TC.PRComp_fit);
            %[TC.effComp_grid_fit] = gridfit(PR,RPM,eff,TC.PRComp_fit,TC.RPMComp_fit);
            Level2 = [0.6:0.05:0.7 0.72:0.02:0.84 0.85];
            Level1 = [0:round(m_dot_max/10/5)*5:m_dot_max];
            figure
            contour(TC.PRComp_fit,TC.RPMComp_fit,TC.m_dotComp_grid_fit',Level1)
            %line(m_dot,PR,RPM,'marker','.','markersize',4,'linestyle','none');
            line(PR,RPM,m_dot,'marker','.','markersize',4,'linestyle','none');

            figure
            %contour(TC.PRComp_fit,TC.RPMComp_fit,TC.effComp_grid_fit,Level2)
            contour(TC.m_dotComp_fit,TC.PRComp_fit,TC.effComp_grid_fit,Level2)
            %line(m_dot,PR,RPM,'marker','.','markersize',4,'linestyle','none');
            line(m_dot,PR,eff,'marker','.','markersize',4,'linestyle','none');
            %line(PR,RPM,eff,'marker','.','markersize',4,'linestyle','none');
            TC.TC_map.comp = struct('RPM',TC.RPMComp_fit','PR',TC.PRComp_fit','eff',TC.effComp_grid_fit, ...
                'm_dot',TC.m_dotComp_grid_fit,'RPM_surge',[0;TC.SAE_map_comp.RPM_surge],...
                'PR_surge',[1;TC.SAE_map_comp.PR_surge],...
                'm_dot_surge',[0;TC.SAE_map_comp.m_dot_surge]);
            prComp_rep = TC.PRComp_fit;
            n288Comp_rep = TC.RPMComp_fit;
            comp_flow_map = TC.m_dotComp_grid_fit;
            comp_eff_map = TC.effComp_grid_fit;
            mDotComp_rep = TC.m_dotComp_fit;
            save(['comp_map' TC.SAE_map_comp.model '.mat'],'prComp_rep','n288Comp_rep','comp_flow_map','comp_eff_map','mDotComp_rep');
            %}

%{
            [xData, yData, zData] = prepareSurfaceData( m_dot, PR, eff );
            % Set up fittype and options.
            ft = 'linearinterp';
            opts = fitoptions( ft );
            opts.Normalize = 'on';
            % Fit model to data.
            TC.fitComp_eff_m_dot_PR = fit( [xData, yData], zData, ft, opts );
            [xData, yData, zData] = prepareSurfaceData( m_dot, PR, RPM);
            TC.fitComp_RPM_m_dot_PR = fit( [xData, yData], zData, ft, opts );
            
            m_dot_min = 0;
            m_dot_max = max(m_dot);
            m_dot_fit = linspace(m_dot_min,m_dot_max,TC.Grid_no+1);
            PR_min = 1.01;
            PR_max = max(PR);
            TC.PRComp_fit = linspace(PR_min,PR_max,TC.Grid_no+1);
            RPM_min = 0;
            RPM_max = max(RPM);
            TC.RPMComp_fit = linspace(RPM_min,RPM_max,TC.Grid_no+1);
            [m_dot_grid, PR_grid] = meshgrid(m_dot_fit,TC.PRComp_fit);
            [RPM_grid, PR_grid] = meshgrid(TC.RPMComp_fit,TC.PRComp_fit);     
            eff_grid = TC.fitComp_eff_m_dot_PR(m_dot_grid, PR_grid);
            RPM_grid_fit = TC.fitComp_RPM_m_dot_PR(m_dot_grid, PR_grid);
                        
            TC.PRComp_fit = TC.PRComp_fit(1:end-1);
            TC.RPMComp_fit = TC.RPMComp_fit(1:end-1);

            TC.effComp_grid_fit = zeros(TC.Grid_no,TC.Grid_no);
            TC.m_dotComp_grid_fit = TC.effComp_grid_fit;
            
            for i=1:TC.Grid_no
                RPM_min = min(RPM_grid_fit(i,:));
                RPM_max = max(RPM_grid_fit(i,:));
                idx_min = find(TC.RPMComp_fit < RPM_min);
                idx_max = find(TC.RPMComp_fit > RPM_max);
                if isempty(idx_max) == 1;
                    rowend = TC.Grid_no;
                else
                    rowend = idx_max(1)-1;
                end;
                if ~isempty(idx_min)
                    for j = idx_min(end)+1:rowend
                        idx_min1 = find(RPM_grid_fit(i,:) < TC.RPMComp_fit(j));
                        idx_min1 = idx_min1(end);
                        if ~isempty(idx_min1)
                            TC.effComp_grid_fit(i,j) = eff_grid(i,idx_min1) + ...
                                (eff_grid(i,idx_min1+1) - eff_grid(i,idx_min1))  / ...
                                (RPM_grid_fit(i,idx_min1+1) - RPM_grid_fit(i,idx_min1)) * ...
                                (TC.RPMComp_fit(j) - RPM_grid_fit(i,idx_min1));
                            TC.m_dotComp_grid_fit(i,j) = m_dot_grid(i,idx_min1) + ...
                                (m_dot_grid(i,idx_min1+1) - m_dot_grid(i,idx_min1))  / ...
                                (RPM_grid_fit(i,idx_min1+1) - RPM_grid_fit(i,idx_min1)) * ...
                                (TC.RPMComp_fit(j) - RPM_grid_fit(i,idx_min1));
                        end;
                    end;
                end;
            end;
            TC.m_dotComp_grid_fit(isnan(TC.m_dotComp_grid_fit)) = 0;
            TC.effComp_grid_fit(isnan(TC.effComp_grid_fit)) = 0;
            % Replace the values at surge area with 0
            %fit_surge = fit([0;TC.SAE_map_comp.RPM_surge],[1;TC.SAE_map_comp.PR_surge],'linearinterp');
            %for i=1:TC.Grid_no
            %    PR_surge_temp = fit_surge(TC.RPMComp_fit(i));
            %    TC.m_dotComp_grid_fit(TC.PRComp_fit > PR_surge_temp,i) = 0;
            %    TC.effComp_grid_fit(TC.PRComp_fit > PR_surge_temp,i) = 0;
            %end;
            
            figure
            contour(RPM_grid_fit,PR_grid,eff_grid,50);
            figure
            contour(RPM_grid_fit,PR_grid,m_dot_grid,50);
            figure
            contour(RPM_grid(1:end-1,1:end-1),PR_grid(1:end-1,1:end-1),TC.effComp_grid_fit,50);
            hold on
            plot([0;TC.SAE_map_comp.RPM_surge],[1;TC.SAE_map_comp.PR_surge],'b');
            hold off
            figure
            contour(RPM_grid(1:end-1,1:end-1),PR_grid(1:end-1,1:end-1),TC.m_dotComp_grid_fit,50);
            hold on
            plot([0;TC.SAE_map_comp.RPM_surge],[1;TC.SAE_map_comp.PR_surge],'b');
            hold off
            TC.TC_map.comp = struct('RPM',TC.RPMComp_fit','PR',TC.PRComp_fit','eff',TC.effComp_grid_fit, ...
                'm_dot',TC.m_dotComp_grid_fit,'RPM_surge',[0;TC.SAE_map_comp.RPM_surge],...
                'PR_surge',[1;TC.SAE_map_comp.PR_surge],...
                'm_dot_surge',[0;TC.SAE_map_comp.m_dot_surge]);
            prComp_rep = TC.PRComp_fit;
            n288Comp_rep = TC.RPMComp_fit;
            comp_flow_map = TC.m_dotComp_grid_fit;
            comp_eff_map = TC.effComp_grid_fit;
            save(['comp_map' TC.SAE_map_comp.model '.mat'],'prComp_rep','n288Comp_rep','comp_flow_map','comp_eff_map');
            %}
            %% Turbine Calculation
            if ~isempty(SAE_map_turb)
                TC.SAE_map_turb = SAE_map_turb;
                if isempty(TC.SAE_map_turb.RPM)
                    TC.ModelType = 'incomplete';
                else
                    TC.ModelType = 'complete';
                end;            
                if strcmp(TC.ModelType,'incomplete')

                    idx = find(TC.effComp_grid_fit == max(max(TC.effComp_grid_fit)));
                    RPM_eff_max = mean(RPM_grid(idx));
                    m_dot_eff_max = mean(TC.m_dotComp_grid_fit(idx));

                    PR = [1;TC.SAE_map_turb.PR];
                    m_dot = [0;TC.SAE_map_turb.m_dot];
                    T_ref = TC.SAE_map_turb.T_ref;
                    eff_max = TC.SAE_map_turb.eta_max;
                    D_wheel = TC.SAE_map_turb.D_wheel;

                    [xData, yData] = prepareCurveData( PR, m_dot );
                    ft = fittype( 'power2' );
                    opts = fitoptions( ft );
                    opts.Display = 'Off';
                    opts.Lower = [-Inf -Inf -Inf];
                    opts.StartPoint = [0.103076329552457 1.24957234727597 -0.0271103634753231];
                    opts.Upper = [Inf Inf Inf];
                    fit_m_dot_PR = fit( xData, yData, ft, opts );

                    % Plot fit with data.
                    figure( 'Name', 'Turbine m_dot vs. PR' );
                    h = plot( fit_m_dot_PR, xData, yData );
                    legend( h, 'm_dot vs. PR', 'fit', 'Location', 'NorthEast' );
                    % Label axes
                    xlabel( 'PR' );
                    ylabel( 'm_dot_corr' );
                    grid on                
        
                    coeff_turb_flow = coeffvalues(fit_m_dot_PR); % a*x^b + c;
                    TC.coeff_m_dot_turb = coeff_turb_flow;
                    save('coeff_turb_flow','coeff_turb_flow');
                    fs = 0.0683;
                    F = 0.4;
                    T_in = 630;
                    Cp_in = 1212.4;
                    K_in = 1.4;
                    m_dot_turb_eff_max = m_dot_eff_max*(1+F*fs);
                    PR_eff_max = 4;
                    err_PR = 1;
                    while err_PR > 0.001
                        m_dot_temp = fit_m_dot_PR(PR_eff_max);
                        f = m_dot_temp - m_dot_turb_eff_max/(PR_eff_max*100)*sqrt(T_in);
                        df = coeff_turb_flow(1)*(1+coeff_turb_flow(2))*PR_eff_max^coeff_turb_flow(2) + coeff_turb_flow(3);
                        dPR = max(-f/df,-0.9*PR_eff_max);
                        err_PR = abs(dPR) / PR_eff_max;
                        PR_eff_max = PR_eff_max + dPR;
                    end;
                    % Fitting eff map
                    TC.Uc_opt = RPM_eff_max*pi/30*D_wheel/2/(sqrt(2*Cp_in*T_in*(1-PR_eff_max^((1-K_in)/K_in))));
                    TC.TC_map.turb = struct('RPM_max',max(RPM),'PR_max',max(PR)','eff_max',eff_max,'Uc_opt',TC.Uc_opt, ...
                        'm_dot_max',max(m_dot),'eta_m',[],'flow_fit_coeff',TC.coeff_m_dot_turb,'D_wheel',D_wheel,'T_ref',T_ref);
                else
                    [~,TC.SAE_map_turb] = TurbMapGenerator(TC.SAE_map_turb,turbRef);
                    m_dot = [0;TC.SAE_map_turb.m_dot];
                    PR = [1;TC.SAE_map_turb.PR];
                    RPM = [0;TC.SAE_map_turb.RPM];
                    eff = [0;TC.SAE_map_turb.eff];
                    m_dot_min = 0;
                    m_dot_max = max(m_dot);
                    TC.m_dotTurb_fit = linspace(m_dot_min,m_dot_max,TC.Grid_no+1);
                    PR_min = 1.01;
                    PR_max = max(PR);
                    TC.PRTurb_fit = linspace(PR_min,PR_max,TC.Grid_no+1);
                    RPM_min = 0;
                    RPM_max = max(RPM);
                    TC.RPMTurb_fit = linspace(RPM_min,RPM_max,TC.Grid_no+1);
                    [TC.m_dotTurb_grid_fit] = gridfit(RPM,PR,m_dot,TC.RPMTurb_fit,TC.PRTurb_fit);
                    [TC.effTurb_grid_fit] = gridfit(RPM,PR,eff,TC.RPMTurb_fit,TC.PRTurb_fit);
                    Level2 = [0.6:0.05:0.7 0.72:0.02:0.84 0.85:0.01 0.9];
                    Level1 = [0:round(m_dot_max/10/5)*5:m_dot_max];
                    figure
                    contour(TC.PRTurb_fit,TC.RPMTurb_fit,TC.m_dotTurb_grid_fit',Level1)
                    %line(m_dot,PR,RPM,'marker','.','markersize',4,'linestyle','none');
                    line(PR,RPM,m_dot,'marker','.','markersize',4,'linestyle','none');

                    figure
                    contour(TC.PRTurb_fit,TC.RPMTurb_fit,TC.effTurb_grid_fit',Level2)
                    %line(m_dot,PR,RPM,'marker','.','markersize',4,'linestyle','none');
                    line(PR,RPM,eff,'marker','.','markersize',4,'linestyle','none');
                    TC.TC_map.turb = struct('RPM',TC.RPMTurb_fit','PR',TC.PRTurb_fit','eff',TC.effTurb_grid_fit, ...
                        'm_dot',TC.m_dotTurb_grid_fit);
                    prTurb_rep = TC.PRTurb_fit;
                    n288Turb_rep = TC.RPMTurb_fit;
                    turb_flow_map = TC.m_dotTurb_grid_fit;
                    turb_eff_map = TC.effTurb_grid_fit;
                    save(['turb_map' TC.SAE_map_comp.model '.mat'],'prTurb_rep','n288Turb_rep','turb_flow_map','turb_eff_map');

                end;
                    
            end;
        end
        function DrawComp(TC,m_dot_user,pr_user,Level)
            RPM = [0;TC.SAE_map_comp.RPM];
            PR = [1;TC.SAE_map_comp.PR];
            m_dot = [0;TC.SAE_map_comp.m_dot];
            m_dot_surge = [0;TC.SAE_map_comp.m_dot_surge];
            m_dot_chock = [0;TC.SAE_map_comp.m_dot_chock];
            PR_surge = [1;TC.SAE_map_comp.PR_surge];
            PR_chock = [1;TC.SAE_map_comp.PR_chock];
            plot(m_dot_surge,PR_surge,'r');
            hold on
            plot(m_dot_chock,PR_chock,'r');
            hold on
            RPM_rep = unique(RPM);
            m = length(RPM_rep);
            for i=1:m
                idx = find(RPM == RPM_rep(i));
                plot(m_dot(idx),PR(idx));
                hold on
            end;
            m_dot_min = 0;
            m_dot_max = max(m_dot);
            m_dot_fit = linspace(m_dot_min,m_dot_max,TC.Grid_no);
            [m_dot_grid, PR_grid] = meshgrid(m_dot_fit,TC.PRTurb_fit);
            eff_grid = TC.fitComp_eff_m_dot_PR(m_dot_grid, PR_grid);
            if isempty(Level)
                Level = [0.6 0.65 0.68 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77];
            end;
            [C,h] = contour(m_dot_grid,PR_grid,eff_grid,Level);
            set(h,'ShowText','on','TextStep',get(h,'LevelStep')/2) 
            hold on
            if ~isempty(m_dot_user)
                plot(m_dot_user,pr_user,'*');
            end;
            hold off
        end
        
        function DrawTurb(TC,pr_user,m_dot_user)
            if  strcmp(TC.ModelType,'incomplete')
                PR = [1;TC.SAE_map_turb.PR];
                m_dot = [0;TC.SAE_map_turb.m_dot];
                plot(PR,m_dot);
                hold on
                plot(pr_user,m_dot_user,'*');
                hold off
            else
                RPM = [0;TC.SAE_map_turb.RPM];
                PR = [1;TC.SAE_map_turb.PR];
                m_dot = [0;TC.SAE_map_turb.m_dot];
                effTurb = [0;TC.SAE_map_turb.eff];
                %m_dot_surge = [0;TC.SAE_map_comp.m_dot_surge];
                %m_dot_chock = [0;TC.SAE_map_comp.m_dot_chock];
                %PR_surge = [1;TC.SAE_map_comp.PR_surge];
                %PR_chock = [1;TC.SAE_map_comp.PR_chock];
                %plot(m_dot_surge,PR_surge,'r');
                %hold on
                %plot(m_dot_chock,PR_chock,'r');
                %hold on
                RPM_rep = unique(RPM);
                m = length(RPM_rep);
                for i=1:m
                    idx = find(RPM == RPM_rep(i));
                    h = plotyy(PR(idx),m_dot(idx),PR(idx),effTurb(idx));
                    xlim(h(1),[1 2.5]);
                    xlim(h(2),[1 2.5]);
                    h(1).XTick = 1:0.2:2.5;
                    h(2).XTick = 1:0.2:2.5;
                    ylim(h(1),[0 2*max(m_dot)]);
                    ylim(h(2),[0 0.8]);
                    %h(1).YTick = 0:0.1:2.5;
                    h(2).YTick = 0:0.1:0.8;
                                        hold on
                end;
            end;
        end

        function Eff = GetEffComp_m_dot_PR(TC,m_dot,PR)
            Eff = TC.fitComp_eff_m_dot_PR(m_dot,PR);
        end
        function Eff = GetEffTurb(TC,ER,T_in,F_in,RPM)
            if TC.ModelType == 'incomplete'
                U1 = 0.5*TC.SAE_map_turb.D_wheel*RPM*pi/30;
                T_in = T_in/1000;
                Cp = (TC.Cp_pTF_coeff * [ER;1])'*([1;T_in;F_in;T_in^2;T_in*F_in;F_in^2]);
                K = (TC.K_pTF_coeff * [ER;1])'*([1;T_in;F_in;T_in^2;T_in*F_in;F_in^2]);
                U2 = sqrt(2*Cp*T_in*1000*(1-ER^((1-K)/K)));
                Uc = U1/U2/TC.Uc_opt;
                Eff = TC.SAE_map_turb.eta_max*(2*Uc-Uc^2);
            else
            end;
        end    
        function RPM = GetRPMComp(TC,m_dot,PR,T_in)
            RPM = TC.fitComp_RPM_m_dot_PR(m_dot,PR).*sqrt(T_in./298.15);
        end
        function m_dot = GetFlowComp_RPM_PR(TC,RPM,PR)
            m_dot = interp2(TC.RPMComp_fit,TC.PR_fit,TC.m_dotComp_grid_fit,RPM,PR);
        end
        function T_out = GetToutComp_RPM_PR(TC,RPM,PR,T_in)
           % Only estimate (+-5%)
            Eff = interp2(TC.RPMComp_fit,TC.PR_fit,TC.effComp_grid_fit,RPM,PR);
            T_out = T_in*(1+(PR^(0.2857)-1)/Eff);
        end
        function m_dot = GetFlowTurb_PR(TC,p_in,p_out,T_in)
        % p_in, p_out in Pa and T_in K
            ER = p_in ./ p_out;
            m_dot_corr = TC.TC_map.turb.flow_fit_coeff(1)*ER.^TC.TC_map.turb.flow_fit_coeff(2) + TC.TC_map.turb.flow_fit_coeff(3);
            m_dot = m_dot_corr.*p_in/1000./sqrt(T_in);
        end
        function ER = GetERTurb_flow(TC, dm_turb, p_out, T_in)
        % dm_turb in kg/s, p_out in Pa and T_in K
            a = TC.TC_map.turb.flow_fit_coeff;
            m = length(dm_turb);
            n = length(p_out);
            if n == 1
                p_out = p_out*ones(m,1);
            end
            p_in = 4*p_out;
            for i = 1:m
                err_dm = 1;
                j = 0;
                while err_dm > 0.001
                    j = j + 1;
                    m_dot_temp = TC.GetFlowTurb_PR(p_in(i),p_out(i),T_in(i));
                    f = m_dot_temp - dm_turb(i);
                    df = (a(1)*(1+a(2))*(p_in(i) / p_out(i))^a(2) + a(3))/sqrt(T_in(i));
                    dp_in = max(-f/df,-0.9*p_in(i));
                    err_dm = abs(f)/dm_turb(i);
                    p_in(i) = p_in(i) + dp_in;
                    if j > 10000
                        break
                        fprintf('Iteration number out of bound (10000) for %ith element',i);
                    end;
                end;            
            end;
            ER = p_in./p_out;
        end
        function WriteTCMap(TC) %Create *.dat file for TC_map
            
            fid1 = fopen('comp_flow_map.dat','w');
            fid2 = fopen('comp_eff_map.dat','w');
            fprintf(fid1,'%u \n',TC.Grid_no);
            fprintf(fid2,'%u \n',TC.Grid_no);
            for i = 1:TC.Grid_no
                fprintf(fid1,'%12.2f',TC.TC_map.comp.RPM(i));
                fprintf(fid2,'%12.2f',TC.TC_map.comp.RPM(i));
            end;
            fprintf(fid1,'\n');
            fprintf(fid2,'\n');            
            for i = 1:TC.Grid_no
                fprintf(fid1,'%12.6f',TC.TC_map.comp.PR(i));
                fprintf(fid2,'%12.6f',TC.TC_map.comp.PR(i));
            end;
            fprintf(fid1,'\n');
            fprintf(fid2,'\n');            
            for i = 1:TC.Grid_no
                for j = 1:TC.Grid_no
                    fprintf(fid1,'%12.6f',TC.TC_map.comp.m_dot(i,j));
                    fprintf(fid2,'%12.6f',TC.TC_map.comp.eff(i,j));
                end
                fprintf(fid1,'\n');
                fprintf(fid2,'\n');                            
            end;            
            fclose(fid1);
            fclose(fid2);
        end        
    end
end
