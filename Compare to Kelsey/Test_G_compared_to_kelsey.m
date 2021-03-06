% From NASA paper: delta_G(T)_species  = G(T)_species - sum(G(T)_elements)
% So for the various reactions, need to read Al(OH)3, AlO(OH), Al2O3, H2O,
% Al, O2, H2. These raw values from 300-600 K are in data folder

%% Reaction Conditions / Constants
T_rxn = 473.15;         % K
T_0 = 298.15;           % K
%P_rxn = 137895.146;
P_rxn = 6.9E+6;         % Pa
P_0 = 101325;           % Pa

R = 8.314;              % J/(mol K)

v_aloh3 = 1/2420;       % m^3/kg
v_alooh = 1/3010;       % m^3/kg
v_al2o3 = 1/3950;       % m^3/kg
v_h2o = 1/1000;         % m^3/kg
v_al = 1/2700;          % m^3/kg

%% For pressure iteration
P_start = P_0;
P_end = 1E+7; % Pa
P_steps = 100;
P = linspace(P_start, P_end, P_steps);

%% Data
addpath('data');
al_raw_data = csvread('al_nasa_raw.csv', 2, 0);
h2_raw_data = csvread('h2_nasa_raw.csv', 2, 0);
o2_raw_data = csvread('o2_nasa_raw.csv', 2, 0);
h2o_raw_data = csvread('h2o_nasa_raw.csv', 2, 0);
aloh3_raw_data = csvread('aloh3_nasa_raw.csv', 2, 0);
% (DEAL WITH THIS LATER) alooh_raw_data = csvread('alooh_nasa_raw.csv', 2, 0);
al2o3_raw_data = csvread('al2o3_nasa_raw.csv', 2, 0);
T = al_raw_data(:, 1);

%% Base elements - G(T)
% Column 6 is -(G-H298)/T, so convert
g_al = -1 * al_raw_data(:, 5) .* al_raw_data(:, 1) + 1E3*(al_raw_data(1, 6) - al_raw_data(1, 3));
g_h2 = -1 * h2_raw_data(:, 5) .* h2_raw_data(:, 1) + 1E3*(h2_raw_data(1, 6) - h2_raw_data(1, 3));
g_o2 = -1 * o2_raw_data(:, 5) .* o2_raw_data(:, 1) + 1E3*(o2_raw_data(1, 6) - o2_raw_data(1, 3));

%% Compounds - G(T) 
g_al2o3 = -1 * al2o3_raw_data(:, 5) .* al2o3_raw_data(:, 1) + 1E3*(al2o3_raw_data(1, 6) - al2o3_raw_data(1, 3));
g_aloh3 = -1 * aloh3_raw_data(:, 5) .* aloh3_raw_data(:, 1) + 1E3*(aloh3_raw_data(1, 6) - aloh3_raw_data(1, 3));
% (DEAL WITH THIS LATER) g_alooh = -1 * alooh_raw_data(:, 5) .* alooh_raw_data(:, 1) + 1E3*(alooh_raw_data(1, 6) - alooh_raw_data(1, 3));
g_h2o = -1 * h2o_raw_data(:, 5) .* h2o_raw_data(:, 1) + 1E3*(h2o_raw_data(1, 6) - h2o_raw_data(1, 3));

%% Compounds - delta_G(T)
delta_g_al2o3 = g_al2o3 - 2*g_al - 3/2*g_o2;
delta_g_aloh3 = g_aloh3 - g_al - 3/2*g_o2 - 3/2*g_h2;
% (DEAL WITH THIS LATER) delta_g_alooh = g_alooh - g_al - g_o2 - 1/2*g_h2;
delta_g_alooh = 1E3.*[-917.916 -904.720 -891.595 -878.351 -865.153 -851.917 -838.740]';
T_alooh = [300 350 400 450 500 550 600]';

%%Not needed in compare to Kelsey Script
%delta_g_h2o = g_h2o - g_h2 - 1/2*g_o2;
% Apply effect of pressure over P range
%delta_g_al2o3 = delta_g_al2o3 + v_al2o3*(P_rxn - P_0);
%delta_g_aloh3 = delta_g_aloh3 + v_aloh3*(P_rxn - P_0);
%delta_g_alooh = delta_g_alooh + v_alooh*(P_rxn - P_0);
%delta_g_h2o = delta_g_h2o + v_h2o*(P_rxn - P_0);
%%

figure(1)
plot(T,delta_g_al2o3);
hold on
plot(T,delta_g_aloh3);
hold on
plot(T_alooh,delta_g_alooh);
legend({'Al_2O_3','Al(OH)_3', 'AlO(OH)'}, 'FontSize', 12);
title('Current Analysis Gibbs Values');

%%
%Plot Kelsey's gibbs values over the same temperature range
x = linspace(300,600,300);
al2o3 = zeros(1,300);
alooh = zeros(1,300);
aloh3 = zeros(1,300);
for n = 1:300
    al2o3(n) = Kelsey_DeltaG_AL2O3_atm(x(n));
    alooh(n) = Kelsey_DeltaG_ALOOH_atm(x(n));
    aloh3(n) = Kelsey_DeltaG_ALOH3_atm(x(n));

    
end
figure(2)
plot (x,al2o3)
hold on
plot (x,aloh3)
hold on 
plot (x,alooh)
legend({'Al_2O_3','Al(OH)_3', 'AlO(OH)'}, 'FontSize', 12);
title('Kelsey Gibbs Values');
%%

