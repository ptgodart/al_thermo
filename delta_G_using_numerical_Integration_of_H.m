%loads matricies of [temp,H] for al2o3 and aloh3
%data ranges from temps 300k-600k in increments of 10k
addpath('data');
H_al2o3_different_temperatures;
H_aloh3_different_temperatures;
H_alooh_different_temperatures;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Al2O3: 
G_al2o3 = zeros(length(H_al2o3),2);
G_al2o3(1) = -1582; %STP value of G_Al2O3 [kj/mol]

%calculate gibbs for Al2O3 at different temps: 
for x  = 2:length(H_al2o3)
    T_a = H_al2o3(x-1,1);
    T_b = H_al2o3(x,1);
    H_a = H_al2o3(x-1,2);
    H_b = H_al2o3(x,2);
    G = T_b*(((T_b- T_a)/6) * ((H_a/(T_a^2)) + 4*(((H_a+ H_b)/2) / (((T_a + T_b)/2)^2)) + (H_b/(T_b^2)))); %Simpson's rule
    G_al2o3(x-1,1) = T_b;
    G_al2o3(x-1,2) = (G_al2o3(x-1)*(T_b/T_a) + G); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Al(OH)3: 

G_aloh3 = zeros(length(H_al2o3),2);
G_aloh3(1) = -1305; %%STP value of G_Al(OH)3 [kj/mol] 

%calculate gibbs for Al(OH)3 at different temps:
for x  = 2:length(H_al2o3)
    T_a = H_aloh3(x-1,1);
    T_b = H_aloh3(x,1);
    H_a = H_aloh3(x-1,2);
    H_b = H_aloh3(x,2);
    G = T_b*(((T_b - T_a)/6) * ((H_a/(T_a^2)) + 4*(((H_a+ H_b)/2) / (((T_a + T_b)/2)^2)) + (H_b/(T_b^2)))); %Simpson's rule
    G_aloh3(x-1,1) = T_b;
    G_aloh3(x-1,2) = (G_aloh3(x-1)*(T_b/T_a) + G);
end 

%AlOOH: 

G_alooh = zeros(length(H_alooh),2);
G_alooh(1) = -918; %%STP value of G_AlOOH [kj/mol] 

%calculate gibbs for AlOOH at different temps:
for x  = 2:length(H_alooh)
    T_a = H_alooh(x-1,1);
    T_b = H_alooh(x,1);
    H_a = H_alooh(x-1,2);
    H_b = H_alooh(x,2);
    G = T_b*(((T_b- T_a)/6) * ((H_a/(T_a^2)) + 4*(((H_a+ H_b)/2) / (((T_a + T_b)/2)^2)) + (H_b/(T_b^2)))); %Simpson's rule
    G_alooh(x-1,1) = T_b;
    G_alooh(x-1,2) = (G_alooh(x-1)*(T_b/T_a) + G); 
end 

G_alooh2 = [300	350	400	450	500	550	600;
    -917.916	-904.72	-891.595	-878.351	-865.153	-851.917	-838.74]';
