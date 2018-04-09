% Reactor specifics
T_rxn = 473.15; % K
T_0 = 298.15;   % K
P_rxn = 6.9E+6; % Pa
P_0 = 101325;   % Pa

% Constants
R = 8.314; % J/(mol K)

delta_g_0_hydroxide = -1151.9E+3;       % J/mol
delta_g_0_oxyhydroxide = -915.04E+3;    % J/mol
delta_g_0_oxide = -1582.4E+3;           % J/mol
delta_g_0_al = 0;
delta_g_0_hydrogen = 0;
delta_g_0_water = -237.14E+3;           % J/mol

delta_h_0_hydroxide = -1277E+3;         % J/mol
delta_h_0_oxyhydroxide = -996E+3;       % J/mol
delta_h_0_oxide = -1669.8E+3;           % J/mol
delta_h_0_al = 0;
delta_h_0_hydrogen = 0;
delta_h_0_water = -285.8E+3;            % J/mol

v_hydroxide = 1/2420;                   % m^3/kg
v_oxyhydroxide = 1/3010;                % m^3/kg
v_oxide = 1/3950;                       % m^3/kg
v_water = 1/1000;                       % m^3/kg
v_al = 1/2700;                          % m^3/kg

delta_G_using_numerical_Integration_of_H;

% Reaction 1: 2Al + 6H2O ==> 2Al(OH)3 + 3H2
% Reaction 2: 2Al + 4H2O ==> 2AlO(OH) + 3H2
% Reaction 3: 2Al + 3H2O ==> Al2O3 + 3H2
% T_array = 300:1:2000;
% delta_G_1_array = zeros(numel(T_array), 2);
% delta_G_2_array = zeros(1, numel(T_array));
% delta_G_3_array = zeros(1, numel(T_array));
% Loop over temperatures in reasonable range, compute G at each temp
% for i = 1:numel(T_array)
%     T = T_array(i);
%     % Gibbs-Helmholtz relation, assuming water, aluminum, and aluminum
%     % byproducts are all incompressible:
%     delta_G_hydroxide = T*(delta_g_0_hydroxide/T_0 + ...
%         delta_h_0_hydroxide*(1/T + 1/T_0)) + v_hydroxide*(P_rxn - P_0);
%     delta_G_oxyhydroxide = T*(delta_g_0_oxyhydroxide/T_0 + ...
%         delta_h_0_oxyhydroxide*(1/T + 1/T_0)) + v_oxyhydroxide*(P_rxn - P_0);
%     delta_G_oxide = T*(delta_g_0_oxide/T_0 + ...
%         delta_h_0_oxide*(1/T + 1/T_0)) + v_oxide*(P_rxn - P_0);

%     delta_G_water = T*(delta_g_0_water/T_0 + ...
%         delta_h_0_water*(1/T + 1/T_0)) + v_water*(P_rxn - P_0);
%     delta_G_al = T*(delta_g_0_al/T_0 + ...
%         delta_h_0_al*(1/T + 1/T_0)) + v_al*(P_rxn - P_0);
%     
%     % Using ideal gas relation for H2:
%     delta_G_h2 = T*(delta_g_0_hydrogen/T_0 + ...
%         delta_h_0_hydrogen*(1/T + 1/T_0)) + R*T*log(P_rxn/P_0);
%     
%     % Reaction 1 (apply stoich ratios)
%     delta_G_products_1 = 2*delta_G_hydroxide + 3*delta_G_h2;
%     delta_G_reactants_1 = 6*delta_G_water + 2*delta_G_al;
%     delta_G_1 = delta_G_products_1 - delta_G_reactants_1;
%     
%     % Reaction 2 (apply stoich ratios)
%     delta_G_products_2 = 2*delta_G_oxyhydroxide + 3*delta_G_h2;
%     delta_G_reactants_2 = 4*delta_G_water + 2*delta_G_al;
%     delta_G_2 = delta_G_products_2 - delta_G_reactants_2;
%  
%     % Reaction 3 (apply stoich ratios)
%     delta_G_products_3 = delta_G_oxide + 3*delta_G_h2;
%     delta_G_reactants_3 = 3*delta_G_water + 2*delta_G_al;
%     delta_G_3 = delta_G_products_3 - delta_G_reactants_3;
%     
%     % Collect values in array
%     delta_G_1_array(i) = delta_G_1;
%     delta_G_2_array(i) = delta_G_2;
%     delta_G_3_array(i) = delta_G_3;
% end

%% Al(OH)_3

for i=1:size(G_aloh3, 1)
    T = G_aloh3(i,1);
    delta_G_hydroxide = G_aloh3(i,2);
    delta_G_water = G_h2o(i,2)
    %T*(delta_g_0_water/T_0 + ...
    %    delta_h_0_water*(1/T + 1/T_0)) + v_water*(P_rxn - P_0);
    delta_G_al = T*(delta_g_0_al/T_0 + ...
        delta_h_0_al*(1/T + 1/T_0)) + v_al*(P_rxn - P_0)
    
    % Using ideal gas relation for H2:
    delta_G_h2 = T*(delta_g_0_hydrogen/T_0 + ...
        delta_h_0_hydrogen*(1/T + 1/T_0)) + R*T*log(P_rxn/P_0);
    
    % Reaction 1 (apply stoich ratios)
    delta_G_products_1 = 2*delta_G_hydroxide + 3*delta_G_h2;
    delta_G_reactants_1 = 6*delta_G_water + 2*delta_G_al;
    delta_G_1 = delta_G_products_1 - delta_G_reactants_1;
    
    delta_G_1_array(i, 1) = T;
    delta_G_1_array(i, 2) = delta_G_1;
end

%% AlOOH

for i=1:size(G_alooh, 1)
    T = G_alooh(i,1);
    delta_G_oxyhydroxide = G_alooh(i,2);
    delta_G_water = G_h2o(i,2)
    %T*(delta_g_0_water/T_0 + ...
    %    delta_h_0_water*(1/T + 1/T_0)) + v_water*(P_rxn - P_0);
    delta_G_al = T*(delta_g_0_al/T_0 + ...
        delta_h_0_al*(1/T + 1/T_0)) + v_al*(P_rxn - P_0);
    
    % Using ideal gas relation for H2:
    delta_G_h2 = T*(delta_g_0_hydrogen/T_0 + ...
        delta_h_0_hydrogen*(1/T + 1/T_0)) + R*T*log(P_rxn/P_0);
    
    % Reaction 2 (apply stoich ratios)
    delta_G_products_2 = 2*delta_G_oxyhydroxide + 3*delta_G_h2;
    delta_G_reactants_2 = 4*delta_G_water + 2*delta_G_al;
    delta_G_2 = delta_G_products_2 - delta_G_reactants_2;
    
    delta_G_2_array(i, 1) = T;
    delta_G_2_array(i, 2) = delta_G_2;
end

%% Al_2O_3

for i=1:size(G_al2o3, 1)
    T = G_al2o3(i,1);
    delta_G_oxide = G_al2o3(i,2);
    delta_G_water = G_h2o(i,2); 
    %T*(delta_g_0_water/T_0 + ...
    %    delta_h_0_water*(1/T + 1/T_0)) + v_water*(P_rxn - P_0);
    delta_G_al = T*(delta_g_0_al/T_0 + ...
        delta_h_0_al*(1/T + 1/T_0)) + v_al*(P_rxn - P_0);
    
    % Using ideal gas relation for H2:
    delta_G_h2 = T*(delta_g_0_hydrogen/T_0 + ...
        delta_h_0_hydrogen*(1/T + 1/T_0)) + R*T*log(P_rxn/P_0);
    
    % Reaction 3 (apply stoich ratios)
    delta_G_products_3 = delta_G_oxide + 3*delta_G_h2;
    delta_G_reactants_3 = 3*delta_G_water + 2*delta_G_al;
    delta_G_3 = delta_G_products_3 - delta_G_reactants_3;
    
    delta_G_3_array(i, 1) = T;
    delta_G_3_array(i, 2) = delta_G_3;
end

% Plotting
figure(1);
clf;
hold on;
plot(delta_G_1_array(:,1), delta_G_1_array(:,2))
plot(delta_G_2_array(:,1), delta_G_2_array(:,2))
plot(delta_G_3_array(:,1), delta_G_3_array(:,2))
legend('Hydroxide', 'Oxyhydroxide', 'Oxide');