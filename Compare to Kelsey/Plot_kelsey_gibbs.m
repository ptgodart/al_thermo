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
