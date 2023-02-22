function [MRho,MQ,jam_speed,mt,mx] = ...
    get_macro_filter_data_ring(times,ring_length,veh_X,veh_V)
%% Extract macroscopic flow property fields
t_begin = times(1);
t_end = times(end);
dt = times(2) - times(1); %Assumes uniform time-step

veh_n = length(veh_X(:,1)); %number of vehicles

veh_X_circ = mod(veh_X,ring_length);

wx = 10; % width of kernel in x
wt = 10; % width of kernel in t

mt = t_begin:5:t_end; %Time samples
mx = linspace(0,ring_length,floor(ring_length)); %position samples

[MT,MX] = meshgrid(mt,mx);

MRho = MT*0;% Mesh for densities
MQ = MT*0;% Mesh for flows

disp('Running macro filter...')
%% 
for j = 1:veh_n % loop over all vehicles trajectories
    percent_complete = 100*(j/veh_n);
    clc
    fprintf(strcat('Percent complete:',num2str(percent_complete)))
%     if floor(j/veh_n*50)>floor((j-1)/veh_n*50) % progress marker
%         fprintf('^')
%     end
    for i = 1:length(times) % loop over all points in that trajectory
        
        dist = min(abs(MX-veh_X_circ(j,i)),abs(MX-veh_X_circ(j,i)-ring_length));
        dist = min(dist,abs(MX-veh_X_circ(j,i)+ring_length));
        
        % Apply a guassian kernel to smooth:
        G = exp(-(dist/wx).^2-...
            ((MT-times(i))/wt).^2)/(pi*wx*wt)*dt; % Gaussian kernel
        
        MRho = MRho+G; % update vehicle density field
        MQ = MQ+G*veh_V(j,i); % update flow rate field
    end
end
fprintf('\n')
disp('finished with kernel filtering.')
%% LSTSQ Rrgression:
fprintf('\n')

jam_speed = zeros(size(mt));

for i=1:length(mt)
    %At each time sample estimates the traveling wave speed.
    jam_line = polyfit(1000*MRho(:,i),3600*MQ(:,i),1);
    jam_speed(i) = jam_line(1);
end
disp('Jam line found.')