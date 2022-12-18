%Benjamin Henry 
%CEE586: Physical Hydrology
%Overland Flow Assigment 
%
%
% This script uses the explicit, finite difference approximation to 
% simulate runoff over an impervious hillslope for specfic storm conditions
% assuming that the kinematic wave approximation is valid 


% Preliminaries 
H_L = 100;          % hillslope length (meters)
S_0 = 0.05;         % Hill Slope (unitless)
S_f = S_0;          % kinematic wave approximation is valid
P_i = 0;            % initial pressure (meters)
M_n = 0.036;        % Manning's number (seconds/m^(1/3))
T_r = 1080;         % Rain Duration (s) = (0.3 h)
R_r = 0.00000556;   % Rain Rate (meters/second)
T_tot = 3600;          % Total Simulation Time (seconds) 

% Moving On
dx = 10;            % delta x (meters)
N = H_L/dx; % Number of Iterations in Space
dt = 1;     % I'm not sure how to use CLT criteria if v is unknown. 
            % I'm just making it 1 second  for now
M = T_tot/dt; % Number of time iterations


%Create an Array for Rain Duration
R = zeros(length(M));

for index = 1:M+1
    if index >= M/2
        R(index) = 0;
    else
         R(index) = R_r;
    end
end
   
%Boundary Conditions/Initial Value 
%Pressure Array Initial Value
for index=1:N
    h(1,index) = 0;
end

% Finite Difference Array
for j=1:M+1
    for i=1:N
        if i < N
        h(j+1,i) = h(j,i) - (dt*sqrt(S_0)/M_n)*...
            (((h(j,i+1).^(5/3))-(h(j,i).^(5/3)))/dx) + R(j)*dt;
        elseif i == N
                h(j+1,i) = h(j,i) - (dt*sqrt(S_0)/M_n)*...
            ((-(h(j,i).^(5/3)))/dx) + R(j)*dt;
        end
    end
end

%Plot Results
figure(1)
plot(h(:,1))
xlabel('Time Iterations')
ylabel('h')
