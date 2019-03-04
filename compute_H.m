%% Data
addpath('data');
h2o_raw_data = csvread('h2o_nasa_raw.csv', 2, 0);
aloh3_raw_data = csvread('aloh3_nasa_raw.csv', 2, 0);
al2o3_raw_data = csvread('al2o3_nasa_raw.csv', 2, 0);
al_raw_data = csvread('al_nasa_raw.csv', 2, 0);
h2_raw_data = csvread('h2_nasa_raw.csv', 2, 0);
o2_raw_data = csvread('o2_nasa_raw.csv', 2, 0);

% H_alooh = [-887.219 -996.389 -996.415 -996.846 -997.269 -997.367 -997.365 -997.058 -996.990]'; % kJ/mol
F_H_alooh = [-4.7715 0 0.335 8.553 15.459 21.375 26.511 31.017 35.001]; % J/mol-K (from Hemingway, 1991)
T_alooh = [273.15 298.15 300 350 400 450 500 550 600]; % K
H_alooh_298 = -996.39E3;
H_alooh = (F_H_alooh .* T_alooh + H_alooh_298)*1E-3;

T_nasa = h2o_raw_data(:, 1);
H_h2o = h2o_raw_data(:, 6);
H_aloh3 = aloh3_raw_data(:, 6);
H_al2o3 = al2o3_raw_data(:, 6);
H_al = al_raw_data(:, 6);
H_h2 = h2_raw_data(:, 6);
H_o2 = o2_raw_data(:, 6);

% Interpolate nice values
T_test = [0 25 50 75 100 125 150 175 200 225 250 275 300]' + 273.15;
H_alooh_interp = interp1(T_alooh, H_alooh, T_test);
H_aloh3_interp = interp1(T_nasa, H_aloh3, T_test);
H_al2o3_interp = interp1(T_nasa, H_al2o3, T_test);
H_h2o_interp = interp1(T_nasa, H_h2o, T_test);
H_al_interp = interp1(T_nasa, H_al, T_test);
H_h2_interp = interp1(T_nasa, H_h2, T_test);
H_o2_interp = interp1(T_nasa, H_o2, T_test);

%% Applying stoichiometry
% Reaction 1: 2Al + 6H2O -> 3H2 + 2Al(OH)3
H_rxn_1 = 3*H_h2_interp + 2*H_aloh3_interp - 6*H_h2o_interp - 2*H_al_interp;
% Reaction 2: 2Al + 4H2O -> 3H2 + 2AlOOH
H_rxn_2 = 3*H_h2_interp + 2*H_alooh_interp - 4*H_h2o_interp - 2*H_al_interp;
% Reaction 3: 2Al + 3H2O -> 3H2 + Al2O3
H_rxn_3 = 3*H_h2_interp + H_al2o3_interp - 3*H_h2o_interp - 2*H_al_interp;

%% Display results
output = [T_test-273.15 H_rxn_1 H_rxn_2 H_rxn_3]; % C kJ/mol kJ/mol kJ/mol
disp(output);