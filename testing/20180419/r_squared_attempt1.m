load('analysis.mat')
 hp_spec = hp_spec-hp_spec(1868);
 x100c_spec = x100c_spec-x100c_spec(1868);
 x50c_spec = x50c_spec-x50c_spec(1868);
 x4c_spec = x4c_spec-x4c_spec(1868);

corrcoef(hp_spec,x100c_spec)
corrcoef(hp_spec,x50c_spec)
corrcoef(hp_spec,x4c_spec)

Rsq1 = 1 - sum((x100c_spec - hp_spec).^2)/sum((x100c_spec - mean(x100c_spec)).^2)
Rsq2 = 1 - sum((x50c_spec - hp_spec).^2)/sum((x50c_spec - mean(x50c_spec)).^2)
Rsq3 = 1 - sum((x4c_spec - hp_spec).^2)/sum((x4c_spec - mean(x4c_spec)).^2)


figure(1)
plot(hp_freq,hp_spec)
hold on
%figure(2)
plot(x100c_freq,x100c_spec)
hold on
%figure(3)
plot(x50c_freq,x50c_spec)
hold on
%figure(4)
plot(x4c_freq,x4c_spec)