% From NASA paper: delta_G(T)_species  = G(T)_species - sum(G(T)_elements)
% So for the various reactions, need to read Al(OH)3, AlO(OH), Al2O3, H2O,
% Al, O2, H2. These raw values from 300-600 K are in data folder

%% Reaction Conditions / Constants
T_rxn = 473.15;         % K
T_0 = 298.15;           % K
%P_rxn = 137895.146;
P_rxn = 0.1E+6;         % Pa
P_0 = 101325;           % Pa

R = 8.314;              % J/(mol K)

v_aloh3 = 1/2420;       % m^3/kg
v_alooh = 1/3010;       % m^3/kg
v_al2o3 = 1/3950;       % m^3/kg
v_h2o = 1/1000;         % m^3/kg
v_al = 1/2700;          % m^3/kg

M_al = 26.9815E-3;      % kg/mol
M_alooh = 59.988E-3;    % kg/mol
M_aloh3 = 78E-3;        % kg/mol
M_al2o3 = 101.96E-3;    % kg/mol
M_h2o = 18.01528E-3;    % kg/mol

%% For pressure iteration
P_start = P_0;
P_end = 1E+7; % Pa
P_steps = 100;
P = linspace(P_start, P_end, P_steps); % Pa

%% Data
addpath('data');
al_raw_data = csvread('al_nasa_raw.csv', 2, 0);
h2_raw_data = csvread('h2_nasa_raw.csv', 2, 0);
o2_raw_data = csvread('o2_nasa_raw.csv', 2, 0);
h2o_raw_data = csvread('h2o_nasa_raw.csv', 2, 0);
h2o_steam_raw_data = csvread('h2o_steam_nasa_raw.csv', 2, 0);
aloh3_raw_data = csvread('aloh3_nasa_raw.csv', 2, 0);
% (DEAL WITH THIS LATER: Issue is that NASA data is only for gas phase HALO2) alooh_raw_data = csvread('alooh_nasa_raw.csv', 2, 0);
al2o3_raw_data = csvread('al2o3_nasa_raw.csv', 2, 0);
T = al_raw_data(:, 1); % K

%% Base elements - G(T)
% Column 6 is -(G-H298)/T, so convert: g [J/mol] = -1 * Value [J/mol-K] * T [K] + 1000*(H [kJ/mol] - (H [kJ/mol] - H298 [kJ/mol])) 
g_al = -1 * al_raw_data(:, 5) .* al_raw_data(:, 1) + 1E3*(al_raw_data(1, 6) - al_raw_data(1, 3)); % J/mol
g_h2 = -1 * h2_raw_data(:, 5) .* h2_raw_data(:, 1) + 1E3*(h2_raw_data(1, 6) - h2_raw_data(1, 3)); % J/mol
g_o2 = -1 * o2_raw_data(:, 5) .* o2_raw_data(:, 1) + 1E3*(o2_raw_data(1, 6) - o2_raw_data(1, 3)); % J/mol

%% Compounds - G(T) 
g_al2o3 = -1 * al2o3_raw_data(:, 5) .* al2o3_raw_data(:, 1) + 1E3*(al2o3_raw_data(1, 6) - al2o3_raw_data(1, 3)); % J/mol
g_aloh3 = -1 * aloh3_raw_data(:, 5) .* aloh3_raw_data(:, 1) + 1E3*(aloh3_raw_data(1, 6) - aloh3_raw_data(1, 3)); % J/mol
g_h2o = -1 * h2o_raw_data(:, 5) .* h2o_raw_data(:, 1) + 1E3*(h2o_raw_data(1, 6) - h2o_raw_data(1, 3)); % J/mol
g_d2o = -1 * h2o_raw_data(:, 5) .* h2o_raw_data(:, 1) + 1E3*(h2o_raw_data(1, 6)*1.03 - h2o_raw_data(1, 3)); % J/mol, multiply H by factor of 1.03 for deuterium
g_h2o_steam = -1 * h2o_steam_raw_data(:, 3) .* h2o_steam_raw_data(:, 1) + 1E3*(h2o_steam_raw_data(1, 4) - h2o_steam_raw_data(1, 2)); % J/mol
%% Compounds - delta_G(T)
delta_g_al2o3 = g_al2o3 - 2*g_al - 3/2*g_o2; % J/mol
delta_g_aloh3 = g_aloh3 - g_al - 3/2*g_o2 - 3/2*g_h2; % J/mol
delta_g_alooh = 1E3.*[-917.916 -904.720 -891.595 -878.351 -865.153 -851.917 -838.740]'; % J/mol, From Hemingway et al, 1991
T_alooh = [300 350 400 450 500 550 600]'; % K
delta_g_h2o = g_h2o - g_h2 - 1/2*g_o2; % J/mol
delta_g_h2o_steam = g_h2o_steam - g_h2 - 1/2*g_o2; % J/mol
% Apply effect of pressure over P range
delta_g_al2o3 = delta_g_al2o3 + v_al2o3*(P - P_0)*M_al2o3; % J/mol
delta_g_aloh3 = delta_g_aloh3 + v_aloh3*(P - P_0)*M_aloh3; % J/mol
delta_g_alooh = delta_g_alooh + v_alooh*(P - P_0)*M_alooh; % J/mol
delta_g_h2o = delta_g_h2o + v_h2o*(P - P_0)*M_h2o; % J/mol (UNCOMMENT WHEN USING NORMAL WATER)
% delta_g_h2o = delta_g_h2o_steam + R*h2o_steam_raw_data(:,1)*log(P/P_0); % ADDED FOR STEAM ANALYSIS. COMMENT FOR NORMAL WATER

%% Elements - delta_G(T)
delta_g_al = zeros(size(T));
delta_g_h2 = zeros(size(T));
% Apply effect of pressure
delta_g_al = delta_g_al + v_al*(P - P_0)*M_al; % J/mol
delta_g_h2 = delta_g_h2 + R*h2_raw_data(:,1)*log(P/P_0); % J/mol

%% Reactions
% REACTION 1: 2Al + 6H2O ==> 2Al(OH)3 + 3H2
delta_g_1 = 2*delta_g_aloh3 + 3*delta_g_h2 - 2*delta_g_al - 6*delta_g_h2o; % J/mol reactant
% REACTION 2: 2Al + 4H2O ==> 2AlO(OH) + 3H2
delta_g_2 = zeros(numel(T_alooh), numel(P));
for i = 1:numel(T_alooh) % Deal with fact that alooh data is sparse
    temp = T_alooh(i);
    for j = 1:numel(T)
        if T(j) == temp
            delta_g_2(i, :) = 2*delta_g_alooh(i, :) + 3*delta_g_h2(j, :) - 2*delta_g_al(j, :) - 4*delta_g_h2o(j, :); % J/mol reactant
            break;
        end
    end
end
% REACTION 3: 2Al + 3H2O ==> Al2O3 + 3H2
delta_g_3 = delta_g_al2o3 + 3*delta_g_h2 - 2*delta_g_al - 3*delta_g_h2o; % J/mol reactant

%% Getting more data
% Fit each delta_g_n with 2nd order polynomial
P_index = round(P_rxn/P_end*P_steps);
T_fit = 273:10:800;
% Reaction 1
delta_g_1_pfit_coeffs = polyfit(T, delta_g_1(:, P_index), 2);
delta_g_1_pfit = polyval(delta_g_1_pfit_coeffs, T_fit);
% Reaction 2
delta_g_2_pfit_coeffs = polyfit(T_alooh, delta_g_2(:, P_index), 2);
delta_g_2_pfit = polyval(delta_g_2_pfit_coeffs, T_fit);
% Reaction 3
delta_g_3_pfit_coeffs = polyfit(T, delta_g_3(:, P_index), 2);
delta_g_3_pfit = polyval(delta_g_3_pfit_coeffs, T_fit);
% Fit surfaces also
[rxn_P_1, rxn_T_1, rxn_delta_g_1] = prepareSurfaceData(P, T, delta_g_1);
[rxn_1_surf_fit, G] = fit([rxn_P_1, rxn_T_1], rxn_delta_g_1, 'cubicinterp');
[P_plot, T_plot] = meshgrid(P, T_fit);

%% Find intersection points
rxn_1_2_intersection = roots(delta_g_1_pfit_coeffs - delta_g_2_pfit_coeffs);
rxn_1_2_intersection = rxn_1_2_intersection(2);
rxn_1_3_intersection = roots(delta_g_1_pfit_coeffs - delta_g_3_pfit_coeffs);
rxn_1_3_intersection = rxn_1_3_intersection(2);
rxn_2_3_intersection = roots(delta_g_2_pfit_coeffs - delta_g_3_pfit_coeffs);
rxn_2_3_intersection = rxn_2_3_intersection(2);

%% Loop through pressures to find transition temps
delta_g_1_2_transitions = zeros(size(P));
delta_g_1_3_transitions = zeros(size(P));
delta_g_2_3_transitions = zeros(size(P));
T_table = [0 25 50 75 100 125 150 175 200 225 250 275 300] + 273.15;
g_table_output = zeros(numel(T_table),numel(P),3);
g_table_1 = zeros(numel(T_table), numel(P));
for p = 1:numel(P)
    % Polyfit
    delta_g_1_coeffs = polyfit(T, delta_g_1(:, p), 2);
    delta_g_2_coeffs = polyfit(T_alooh, delta_g_2(:, p), 2);
    delta_g_3_coeffs = polyfit(T, delta_g_3(:, p), 2);
    % Find intersections
    delta_g_1_2_intersection = roots(delta_g_1_coeffs - delta_g_2_coeffs);
    delta_g_1_2_transitions(p) = delta_g_1_2_intersection(2);
    delta_g_1_3_intersection = roots(delta_g_1_coeffs - delta_g_3_coeffs);
    delta_g_1_3_transitions(p) = delta_g_1_3_intersection(2);
    delta_g_2_3_intersection = roots(delta_g_2_coeffs - delta_g_3_coeffs);
    delta_g_2_3_transitions(p) = delta_g_2_3_intersection(2);
    % Fill in output table for paper
    g_table_output(:, p, 1) = polyval(delta_g_1_coeffs, T_table)';
    g_table_output(:, p, 2) = polyval(delta_g_2_coeffs, T_table)';
    g_table_output(:, p, 3) = polyval(delta_g_3_coeffs, T_table)';
end

%% For paper table
disp(P(1:10:end))
disp([T_table'-273.15 g_table_output(:, 1:10:end, 1)*1E-3]);
disp(P(1:10:end))
disp([T_table'-273.15 g_table_output(:, 1:10:end, 2)*1E-3]);
disp(P(1:10:end))
disp([T_table'-273.15 g_table_output(:, 1:10:end, 3)*1E-3]);

%% Plotting
figure(1); clf;
hold on;
grid on;
T = al_raw_data(:, 1);
plot(T-273.15, delta_g_1(:, P_index));
plot(T_alooh-273.15, delta_g_2(:, P_index));
plot(T-273.15, delta_g_3(:, P_index));
legend({'Al(OH)_3', 'AlO(OH)', 'Al_2O_3'}, 'FontSize', 12);
xlabel('Temperature [ºC]', 'FontSize', 14);
ylabel('Gibbs Free Energy [J/mol-ºC]', 'FontSize', 14);
title(['Gibbs Free Energy For Al-Water Reactions at ' num2str(P(P_index)) ' Pa'], 'FontSize', 16);
% Plot fitted as well
figure(2); clf;
hold on;
grid on;
plot(T_fit-273.15, delta_g_1_pfit);
plot(T_fit-273.15, delta_g_2_pfit);
plot(T_fit-273.15, delta_g_3_pfit);
xlim([0 400]);
legend({'Al(OH)_3', 'AlO(OH)', 'Al_2O_3'}, 'FontSize', 12);
xlabel('Temperature [ºC]', 'FontSize', 14);
ylabel('Extrapolated Gibbs Free Energy [J/mol-ºC]', 'FontSize', 14);
title(['Extrapolated Gibbs Free Energy For Al-Water Reactions at ' num2str(P(P_index)) ' Pa'], 'FontSize', 16);
% Plot surface
figure(3); clf;
axes1 = axes;
hold(axes1,'on');
surf_1 = surf(P/1E6,T-273.15,delta_g_1/1E3, 'DisplayName','Al(OH)_3',...
    'FaceColor',[0 0.447058826684952 0.74117648601532]);
surf_2 = surf(P/1E6,T_alooh-273.15,delta_g_2/1E3, 'DisplayName','AlOOH','FaceColor',[1 0 0]);
surf_3 = surf(P/1E6,T-273.15,delta_g_3/1E3, 'DisplayName','Al_2O_3','FaceColor',[1 1 0]);
xlabel('Pressure [MPa]', 'FontSize', 14);
ylabel('Temperature [ºC]', 'FontSize', 14);
zlabel('Gibbs Free Energy [kJ/mol]', 'FontSize', 14);
title('Gibbs Free Energy For Al-Water Reactions Over Operating Range');
view(axes1,[64 26.7999999999999]);
box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',16,'Location','northeast');
% Plot transitions
figure(4); clf;
hold on;
plot(P, delta_g_1_2_transitions-273.15);
plot(P, delta_g_2_3_transitions-273.15);
%plot(P, delta_g_1_3_transitions-273.15);
plot(P, 1730.63 ./ (8.07131 - log10(0.0075*P)) - 233.426, '--');
title('Predicting Aluminum-Water Reaction Byproducts', 'FontSize', 16);
ylabel('Temperature [ºC]', 'FontSize', 14);
xlabel('Pressure [Pa]', 'FontSize', 14);
xlim([P_start, P_end]);
legend({'Al(OH)_3 --> AlO(OH)', 'AlO(OH) --> Al_2O_3', 'BP of Water'}, 'FontSize', 12);
%legend({'Al(OH)_3 --> AlO(OH)', 'AlO(OH) --> Al_2O_3', 'Al(OH)_3 --> Al_2O_3', 'BP of Water'}, 'FontSize', 12);
