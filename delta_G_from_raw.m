addpath('data');
al_raw_path = "al_nasa_raw.csv";
alooh_raw_path = "alooh_nasa_raw.csv";
aloh3_raw_path = "aloh3_nasa_raw.csv";
al2o3_raw_path = "al2o3_nasa_raw.csv";
h2_raw_path = "h2_nasa_raw.csv";
h2o_raw_path = "h2o_nasa_raw.csv";

al_raw_data = csvread(al_raw_path, 2, 0);
T_al = al_raw_data(:, 1);
G_al = al_raw_data(:, 5);

alooh_raw_data = csvread(alooh_raw_path, 2, 0);
T_alooh = alooh_raw_data(:, 1);
G_alooh = alooh_raw_data(:, 5);

aloh3_raw_data = csvread(aloh3_raw_path, 2, 0);
T_aloh3 = aloh3_raw_data(:, 1);
G_aloh3 = aloh3_raw_data(:, 5);

al2o3_raw_data = csvread(al2o3_raw_path, 2, 0);
T_al2o3 = al2o3_raw_data(:, 1);
G_al2o3 = al2o3_raw_data(:, 5);

h2_raw_data = csvread(h2_raw_path, 2, 0);
T_h2 = h2_raw_data(:, 1);
G_h2 = h2_raw_data(:, 5);

h2o_raw_data = csvread(h2o_raw_path, 2, 0);
T_h2o = h2o_raw_data(:, 1);
G_h2o = h2o_raw_data(:, 5);

range = 60;

%% Al2O3 Reaction
% 2Al + 3H2O ==> Al2O3 + 3H2

for t = 1:range
    temp = T_al2o3(t);
    T_al2o3_rxn(t) = T_al2o3(t);
    G_al2o3_rxn(t) = temp * (3*G_h2(t) + G_al2o3(t) - 2*G_al(t) - 3*G_h2o(t));
end

%% AlO(OH) Reaction
% 2Al + 4H2O ==> 2AlO(OH) + 3H2

for t = 1:range
    temp = T_alooh(t);
    T_alooh_rxn(t) = T_alooh(t);
    G_alooh_rxn(t) = temp * (3*G_h2(t) + 2*G_alooh(t) - 2*G_al(t) - 4*G_h2o(t));
end

%% Al(OH)3 Reaction
% 2Al + 6H2O ==> 2Al(OH)3 + 3H2

for t = 1:range
    temp = T_aloh3(t);
    T_aloh3_rxn(t) = T_aloh3(t);
    G_aloh3_rxn(t) = temp * (3*G_h2(t) + 2*G_aloh3(t) - 2*G_al(t) - 6*G_h2o(t));
end

%% Plotting

figure(1);
clf;
hold on;
plot(T_aloh3_rxn, G_aloh3_rxn);
plot(T_alooh_rxn, G_alooh_rxn);
plot(T_al2o3_rxn, G_al2o3_rxn);
legend('Hydroxide', 'Oxyhydroxide', 'Oxide');