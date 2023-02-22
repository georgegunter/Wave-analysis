%% Define relevant model funtions:

% a=p(1);b=p(2);V0=p(3);delta=p(4);T=p(5);s0=p(6);

IDM_Accel = @(p,s,ds,v) p(1).*( 1 - (v/p(3)).^p(4) - ...
    ((p(6) + max(p(5)*v + v.*(-ds)./(2*sqrt(p(1).*p(2))),0))./s).^2 );

S_eq = @(p,v) (p(6)+p(5).*v)./sqrt(1-(v./p(3)).^p(4));
%% Define the parameters and equilibria:
% of the form: params = [a,b,V0,delta,T,s0]
params = [1.3,1.75,27.0,4.0,2.5,5.0];

% params = [3.6652,0.1550,13.9884,8.3875,-0.0407,8.8713];

v_eq = 1.0:.1:29.5;
s_eq = S_eq(params,v_eq);


%% Plot results:
figure()
plot(s_eq,v_eq)

[q,rho,stab] = string_stability_general(params,IDM_Accel,s_eq,v_eq,4);
set(gca,'FontSize',24)
% close all
%% Plotting:

figure()
subplot(3,1,1)
hold on
hl = plot(v_eq(stab<0),s_eq(stab<0),'r.',v_eq(stab>0),s_eq(stab>0),'b.','MarkerSize',15);
hold off
legend(hl,'unstable','stable','Location','SouthEast')
ylabel('equilibrium gap [m]')
grid on
box on
set(gca,'FontSize',24,'LineWidth',3)
title('IDM equilibrium operating conditions','Fontsize',35)

subplot(3,1,2)
hl = plot(v_eq(stab<0),q(stab<0)*3600,'r.',v_eq(stab>0),q(stab>0)*3600,'b.','MarkerSize',15);
ylabel('flow rate  [veh/hr/lane]')
grid on
set(gca,'FontSize',24,'LineWidth',3)
box on

subplot(3,1,3)
hold on
plot([v_eq(1),v_eq(end)],[0,0],'k--','LineWidth',3)
hl = plot(v_eq(stab<0),stab(stab<0),'r.',v_eq(stab>0),stab(stab>0),'b.','MarkerSize',15);
hold off
xlabel('equilibrium speed / m/s'), ylabel('\lambda (stability value)')
grid on
set(gca,'FontSize',24,'LineWidth',3)
box on