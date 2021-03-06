addpath('data');
al_raw_path = "al_nasa_raw.csv";
alooh_raw_path = "alooh_nasa_raw.csv";
aloh3_raw_path = "aloh3_nasa_raw.csv";
al2o3_raw_path = "al2o3_nasa_raw.csv";
h2_raw_path = "h2_nasa_raw.csv";
h2o_raw_path = "h2o_nasa_raw.csv";

% Table gives values as -(G-H298)/T, so multiply through by T, multiply
% through by -1, and add H298 (which is column 6 minus column 3). Units are
% listed in table.
al_raw_data = csvread(al_raw_path, 2, 0);
H298_al = 1000 * (al_raw_data(1, 6) - al_raw_data(1, 3)); % J/mol
T_al = al_raw_data(:, 1); % K
G_al = (al_raw_data(:, 5) * -1 .* T_al) + H298_al; % J/mol

alooh_raw_data = csvread(alooh_raw_path, 2, 0);
H298_alooh = 1000 * (alooh_raw_data(1, 6) - alooh_raw_data(1, 3));
T_alooh = alooh_raw_data(:, 1);
G_alooh = (alooh_raw_data(:, 5) * -1 .* T_alooh) + H298_alooh;

aloh3_raw_data = csvread(aloh3_raw_path, 2, 0);
H298_aloh3 = 1000 * (aloh3_raw_data(1, 6) - aloh3_raw_data(1, 3));
T_aloh3 = aloh3_raw_data(:, 1);
G_aloh3 = (aloh3_raw_data(:, 5) * -1 .* T_aloh3) + H298_aloh3;

al2o3_raw_data = csvread(al2o3_raw_path, 2, 0);
H298_al2o3 = 1000 * (al2o3_raw_data(1, 6) - al2o3_raw_data(1, 3));
T_al2o3 = al2o3_raw_data(:, 1);
G_al2o3 = (al2o3_raw_data(:, 5) * -1 .* T_al2o3) + H298_al2o3;

h2_raw_data = csvread(h2_raw_path, 2, 0);
H298_h2 = 1000 * (h2_raw_data(1, 6) - h2_raw_data(1, 3));
T_h2 = h2_raw_data(:, 1);
G_h2 = (h2_raw_data(:, 5) * -1 .* T_h2) + H298_h2;

h2o_raw_data = csvread(h2o_raw_path, 2, 0);
H298_h2o = 1000 * (h2o_raw_data(1, 6) - h2o_raw_data(1, 3));
T_h2o = h2o_raw_data(:, 1);
G_h2o = (h2o_raw_data(:, 5) * -1 .* T_h2o) + H298_h2o;

range = 30;

%% Al2O3 Reaction
% 2Al + 3H2O ==> Al2O3 + 3H2
T_al2o3_rxn = zeros(range, 1);
G_al2o3_rxn = zeros(range, 1);
for t = 1:range
    T_al2o3_rxn(t) = T_al2o3(t);
    G_al2o3_rxn(t) = 3*G_h2(t) + G_al2o3(t) - 2*G_al(t) - 3*G_h2o(t);
end

%% AlO(OH) Reaction
% 2Al + 4H2O ==> 2AlO(OH) + 3H2
% Note: this is from NASA thermo builder but only for gaseous alooh
T_alooh_rxn = zeros(range, 1);
G_alooh_rxn = zeros(range, 1);
for t = 1:range
    T_alooh_rxn(t) = T_alooh(t);
    G_alooh_rxn(t) = 3*G_h2(t) + 2*G_alooh(t) - 2*G_al(t) - 4*G_h2o(t);
end

%% Al(OH)3 Reaction
% 2Al + 6H2O ==> 2Al(OH)3 + 3H2
T_aloh3_rxn = zeros(range, 1);
G_aloh3_rxn = zeros(range, 1);
for t = 1:range
    T_aloh3_rxn(t) = T_aloh3(t);
    G_aloh3_rxn(t) = 3*G_h2(t) + 2*G_aloh3(t) - 2*G_al(t) - 6*G_h2o(t);
end

%% AlO(OH) from Hemingway paper
% 2Al + 4H2O ==> 2AlO(OH) + 3H2
H_298 = -996.38E3;
G_alooh_paper_eq = [37.191 37.880 39.485 41.656 44.180 46.922 49.796]; % J/mol-K
% Values right from paper
G_alooh_paper_raw = 1E3*[-917.916 -904.720 -891.595 -878.351 -865.153 -851.917 -838.740];
T_alooh_paper = [300 350 400 450 500 550 600];
% Converted values
G_alooh_paper = -1*T_alooh_paper.*G_alooh_paper_eq + H_298;
T_alooh_paper_rxn = zeros(7,1);
G_alooh_paper_rxn_1 = zeros(7,1);
G_alooh_paper_rxn_2 = zeros(7,1);
for t = 1:7
    T_alooh_paper_rxn(t) = T_alooh_paper(t);
    G_alooh_paper_rxn_1(t) = 3*G_h2(t) + 2*G_alooh_paper(t) - 2*G_al(t) - 4*G_h2o(t);
    G_alooh_paper_rxn_2(t) = 3*G_h2(t) + 2*G_alooh_paper_raw(t) - 2*G_al(t) - 4*G_h2o(t);
end

%% Plotting
figure(1);
clf;
hold on;
plot(T_aloh3_rxn - 274.15, G_aloh3_rxn);
%plot(T_alooh_rxn - 274.15, G_alooh_rxn);
plot(T_al2o3_rxn - 274.15, G_al2o3_rxn);
plot(T_alooh_paper_rxn - 274.15, G_alooh_paper_rxn_1);
%plot(T_alooh_paper_rxn - 274.15, G_alooh_paper_rxn_2);
legend('Hydroxide', 'Oxide', 'Oxyhydroxide - paper (converted)', 'Oxyhydroxide - paper (from table)');