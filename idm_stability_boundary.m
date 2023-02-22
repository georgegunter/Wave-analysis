%% Define relevant model funtions:

% a=p(1);b=p(2);V0=p(3);delta=p(4);T=p(5);s0=p(6);

IDM_Accel = @(p,s,ds,v) p(1).*( 1 - (v/p(3)).^p(4) - ...
    ((p(6) + max(p(5)*v + v.*(-ds)./(2*sqrt(p(1).*p(2))),0))./s).^2 );

S_eq = @(p,v) (p(6)+p(5).*v)./sqrt(1-(v./p(3)).^p(4));
%% Analyze at a certain speed:

v_eq = 5.0;
v0=30;T=1;delta=4;s0=2; %parameters that will stay fixed.
s_eq = S_eq([1.0,1.0,v0,delta,T,s0],v_eq);
    
a_range = .1:.05:3.0;
b_range = a_range;

stab_vals = zeros(numel(a_range),numel(b_range));

for i=1:length(a_range)
    for j=1:length(b_range)
        a = a_range(i);
        b = b_range(j);
        
        params = [a,b,v0,delta,T,s0];
        [q,rho,stab] = string_stability_general(params,IDM_Accel,s_eq,v_eq,0,false);
        
        stab_vals(i,j) = stab < 0;
        
    end
end
        
disp('Finished')
%% Plotting:


figure()
contourf(a_range,b_range,stab_vals)
ylabel('b')
xlabel('a')
set(gca,'FontSize',20)
title('Stability Boundary')








