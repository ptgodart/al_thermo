% From NASA paper: delta_G(T)_species  = G(T)_species - sum(G(T)_elements)
% So for the various reactions, need to read Al(OH)3, AlO(OH), Al2O3, H2O,
% Al, O2, H2. These raw values from 300-600 K are in data folder

%% Data
addpath('data');
al_raw_data = csvread('al_nasa_raw.csv', 2, 0);
h2_raw_data = csvread('h2_nasa_raw.csv', 2, 0);
o2_raw_data = csvread('o2_nasa_raw.csv', 2, 0);
h2o_raw_data = csvread('h2o_nasa_raw.csv', 2, 0);
aloh3_raw_data = csvread('aloh3_nasa_raw.csv', 2, 0);
% (DEAL WITH THIS LATER) alooh_raw_data = csvread('alooh_nasa_raw.csv', 2, 0);
al2o3_raw_data = csvread('al2o3_nasa_raw.csv', 2, 0);

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
delta_g_h2o = g_h2o - g_h2 - 1/2*g_o2;

%% Elements - delta_G(T)
delta_g_al = zeros(size(T));
delta_g_h2 = zeros(size(T));

%% Reactions
% REACTION 1: 2Al + 6H2O ==> 2Al(OH)3 + 3H2
delta_g_1 = 2*delta_g_aloh3 + 3*delta_g_h2 - 2*delta_g_al - 6*delta_g_h2o;
% REACTION 2: 2Al + 4H2O ==> 2AlO(OH) + 3H2
delta_g_2 = zeros(size(T_alooh));
for i = 1:numel(T_alooh) % Deal with fact that alooh data is sparse
    temp = T_alooh(i);
    for j = 1:numel(T)
        if T(j) == temp
            delta_g_2(i, 1) = 2*delta_g_alooh(i) + 3*delta_g_h2(j) - 2*delta_g_al(j) - 4*delta_g_h2o(j);
            break;
        end
    end
end
% REACTION 3: 2Al + 3H2O ==> Al2O3 + 3H2
delta_g_3 = delta_g_al2o3 + 3*delta_g_h2 - 2*delta_g_al - 3*delta_g_h2o;

%% Plotting
figure;
clf;
hold on;
T = al_raw_data(:, 1);
plot(T, delta_g_1);
plot(T_alooh, delta_g_2);
plot(T, delta_g_3);
legend('Al(OH)_3', 'AlO(OH)', 'Al_2O_3');