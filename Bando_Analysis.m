%% Define relevant model funtions:
V = @(p,s) p(3)*(tanh(s./p(4)-p(5))+tanh(p(5)))/(1+tanh(p(5)));
Bando_Accel = @(p,s,ds,v)  p(1)*(p(3)*(tanh(s./p(4)-p(5))+tanh(p(5)))/(1+tanh(p(5)))-v) + p(2)*((ds)./(s.^2));


%% Define the parameters and equilibria:
% of the form: params = [alpha,beta,v0,d0]
% params = [0.6660,21.5975,8.9368,2.2146,2.8150];
% params = [0.7,25,8.9368,2.2146,2.8150];
% params = [0.8,20.0,15.0,1.0,2.0];

params = [.5,20,30.0,12,2];
s_eq = 0:.1:20;
v_eq = V(params,s_eq);
figure()
plot(s_eq,v_eq,'LineWidth',5)
hold on


len = 4.5;
wantplots = true;
%% Plot results:
close all;
clc;

[q,rho,stab] = string_stability_general(params,Bando_Accel,s_eq,v_eq,0.0);
set(gca,'FontSize',24)
