function [q,rho,stab] = string_stability_general(params,accel_func,s_eq,v_eq,len)
% AUTHORS: Benjamin Seibold and George Gunter (3/11/2020)

    acceleration = @(s,ds,v) accel_func(params,s,ds,v);

    v = v_eq; % sample velocity range
    s = s_eq; % corresponding equilibrium ga

    % Calculate additional fundamental diagram variables
    rho = 1./(s+len); % vehicle density
    q = rho.*v; % flow rate

    % Calculate partial derivatives of acceleration function at equilibrium
    ep = 1e-6;
    dads = (acceleration(s+ep,0,v)-acceleration(s-ep,0,v))/(2*ep);
    dadds = (acceleration(s,ep,v)-acceleration(s,-ep,v))/(2*ep);
    dadv = (acceleration(s,0,v+ep)-acceleration(s,0,v-ep))/(2*ep);

    % Calculate linear system coefficients
    alpha1 = dads;
    alpha2 = dadds-dadv;
    alpha3 = dadds;

    % Calculate stability value
    stab = alpha2.^2-alpha3.^2-2*alpha1; % state stable exactlt when positive

    % Calculate growth factors
    w = 1i*[0,10.^linspace(-3,2,200)]'; % values only imaginary axis
    e = w*0+1;
    F = (e*alpha1+w*alpha3)./... % transfer function
        (e*alpha1+w*alpha2+w.^2*(alpha1*0+1));
    growth = max(abs(F));
    growth(stab>0) = nan; % remove stable states

    % Plotting
    clf
    subplot(2,2,1)
    plot(s,v,'k-','LineWidth',2)
    hold on
    hl = plot(s(stab<0),v(stab<0),'r.',s(stab>0),v(stab>0),'b.');
    hold off
    axis([0,max(s),0,max(v)])
    legend(hl,'unstable','stable','Location','SouthEast')
    xlabel('equilibrium gap / m'), ylabel('equilibrium speed / m/s')
    title('Equilbrium: speed vs. gap')
    grid on
    set(gca,'FontSize',15)

    subplot(2,2,2)
    plot([rho,0],[q,0],'k-','LineWidth',2)
    hold on
    hl = plot(rho(stab<0),q(stab<0)*3600,'r.',rho(stab>0),q(stab>0)*3600,'b.');
    hold off
    xlim([0,max(rho)])
    legend(hl,'unstable','stable','Location','NorthEast')
    xlabel('vehicle density / 1/m'), ylabel('flow rate / 1/hr')
    title('Equilbrium: fundamental diagram (flow vs. density)')
    grid on
    set(gca,'FontSize',15)

    subplot(2,2,3)
    plot([0,max(s)],[0,0],'k-')
    hold on
    hl = [plot(s,stab,s,growth,'LineWidth',2);...
        plot(s,alpha1,s,alpha2,s,alpha3,'LineWidth',1)];
    hold off
    xlim([0,max(s)])
    legend(hl,'stability value (stable: >0)','growth factor (veh. to veh.)',...
        '\alpha_1','\alpha_2','\alpha_3')
    xlabel('equilibrium gap / m'), ylabel('stability parameters')
    title('Stability parameters (gap view)')
    grid on
    set(gca,'FontSize',15)

    subplot(2,2,4)
    plot([0,max(rho)],[0,0],'k-')
    hold on
    hl = plot(rho,stab,rho,(growth-1)*100,'LineWidth',2);
    hold off
    xlim([0,max(rho)])
    legend(hl,'stability value','growth factor (% increase)',...
        'Location','NorthWest')
    xlabel('vehicle density / 1/m'), ylabel('stability parameters')
    title('Stability parameters (density view)')
    grid on
    set(gca,'FontSize',15)
    
    [q_crit,q_crit_index] = max(q*3600);
    v_crit = v(q_crit_index);
    
    
    figure()
    hl = plot(q(stab<0)*3600,v(stab<0),'r.',q(stab>0)*3600,v(stab>0),'b.',...
        [0,q_crit],[v_crit,v_crit],'k--','MarkerSize',15);
    hold off
    xlim([0,max(q*3600)]+100)
    legend(hl,'unstable','stable','Critical Speed','Location','NorthEast')
    xlabel('flow rate / 1/hr'), ylabel('vehicle speed m/s')
    title_string = cell2mat(strcat({'Equilbrium: Flow to Speed,  Q_{critical} = '},{num2str(q_crit)},...
        {' V_{critical} = '},{num2str(v_crit)}));
    title(title_string)
    grid on


    figure()
    hl = plot(rho(stab<0)*1000,q(stab<0)*3600,'r.',rho(stab>0)*1000,q(stab>0)*3600,'b .','MarkerSize',15);
    hold off
    legend(hl,'unstable','stable','Location','NorthEast')
    ylabel('flow rate  [veh/hr]'), xlabel('vehicle density [veh/km]')
    title_string = 'Fundamental Diagram';
    title(title_string)
    grid on
    

end

