% Rishav (2020/12/26)

start_time = 0;
stop_time = 5;
dt = 0.1;
time = start_time:dt:stop_time;

% [del q, beta]
state = zeros(6,time);

% Plot setting
set(groot,'defaulttextinterprete','latex');  
set(groot, 'defaultAxesTickLabelInterprete','latex');  
set(groot, 'defaultLegendInterprete','latex');

