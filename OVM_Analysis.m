%% Define relevant model funtions:
V = @(p,s) p(3)/2*(1 - cos(pi*(s-p(4))/(p(5) - p(4))));
OVM_Accel = @(p,s,ds,v)  p(1)*(p(3)/2*(1 - cos(pi*(s-p(4))/(p(5) - p(4))))-v) + p(2)*(ds);


%% Define the parameters and equilibria:
% of the form: params = [alpha,beta,v0,d0]
params = [.1,5,30,2,15];
s_eq = params(4):.1:params(5);
v_eq = V(params,s_eq);

%% Plot results:
figure()
plot(s_eq,v_eq)

[q,rho,stab] = string_stability_general(params,Bando_Accel,s_eq,v_eq,4);
