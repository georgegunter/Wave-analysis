%% Define relevant model funtions:

% a=p(1);b=p(2);V0=p(3);delta=p(4);T=p(5);s0=p(6);

IDM_Accel = @(p,s,ds,v) p(1).*( 1 - (v/p(3)).^p(4) - ...
    ((p(6) + max(p(5)*v + v.*(-ds)./(2*sqrt(p(1).*p(2))),0))./s).^2 );

S_eq = @(p,v) (p(6)+p(5).*v)./sqrt(1-(v./p(3)).^p(4));
%% Define the parameters and equilibria:
% of the form: params = [a,b,V0,delta,T,s0]
param_set_1 = [0.75,2.0,30.0,4.0,1.0,2.0];
param_set_2 = [1.5,1.0,30.0,4.0,1.0,2.0];

len_1 = 5.0;
len_2 = 15.0;

v_eq = 0:.1:30; 
s_eq = S_eq(param_set_1,v_eq);

[q_1,rho_1,stab_1] = string_stability_general(param_set_1,IDM_Accel,s_eq,v_eq,len_1);
[q_2,rho_2,stab_2] = string_stability_general(param_set_2,IDM_Accel,s_eq,v_eq,len_1);


close all

%%
figure()

subplot(1,2,1)
hl = plot(rho_1(stab_1<0),q_1(stab_1<0)*3600,'r-.',rho_1(stab_1>0),q_1(stab_1>0)*3600,'b-.','LineWidth',5);
hold off
legend(hl,'unstable','stable','Location','NorthEast')
ylabel('flow rate  [veh/hr]'), xlabel('vehicle density [veh/m]')
title('a = 0.75, b = 2.0')
grid on
set(gca,'LineWidth',4)
set(gca,'FontSize',24)

subplot(1,2,2)
hl = plot(rho_2(stab_2<0),q_2(stab_2<0)*3600,'r-.',rho_2(stab_2>0),q_2(stab_2>0)*3600,'b-.','LineWidth',5);
hold off
% legend(hl,'unstable','stable','Location','NorthEast')
ylabel('flow rate  [veh/hr]'), xlabel('vehicle density [veh/m]')
title('a = 1.5, b = 2.0')
grid on
set(gca,'LineWidth',4)
set(gca,'FontSize',24)