load('analysis.mat')
peak_hp = zeros(length(hp_spec),1);
peak_100c = zeros(length(hp_spec),1);
peak_50c = zeros(length(hp_spec),1);
peak_4c = zeros(length(hp_spec),1);

min_dist = 75;
min_height = 0.0001;

[pks1,locs1] = findpeaks(smooth(hp_spec,11),hp_freq, 'MinPeakDistance', min_dist, 'Threshold', min_height);
[pks2,locs2] = findpeaks(smooth(x100c_spec,11),x100c_freq, 'MinPeakDistance', min_dist, 'Threshold', min_height);
[pks3,locs3] = findpeaks(smooth(x50c_spec,11),x50c_freq, 'MinPeakDistance', min_dist, 'Threshold', min_height);
[pks4,locs4] = findpeaks(smooth(x4c_spec,11),x4c_freq, 'MinPeakDistance', min_dist, 'Threshold', min_height);

peak_hp(round((locs1'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
peak_100c(round((locs2'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
peak_50c(round((locs3'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
peak_4c(round((locs4'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;

xcorr_hp_100c = max(xcorr(peak_hp,peak_100c,'coeff'));
xcorr_hp_50c = max(xcorr(peak_hp,peak_50c,'coeff'));
xcorr_hp_4c = max(xcorr(peak_hp,peak_4c,'coeff'));

xcorr_4c_100c = max(xcorr(peak_4c,peak_100c,'coeff'));
xcorr_4c_50c = max(xcorr(peak_4c,peak_50c,'coeff'));
xcorr_4c_hp = max(xcorr(peak_4c,peak_hp,'coeff'));

figure(1); clf;
plot(xcorr(peak_hp, peak_50c));

figure(2); clf;
hold on;
scatter(hp_freq, peak_100c);
scatter(hp_freq, peak_hp);
plot(x100c_freq, x100c_spec);
plot(x100c_freq, hp_spec);
legend('100C peaks', 'HP peaks');

figure(3); clf;
hold on;
scatter(hp_freq, peak_4c);
scatter(hp_freq, peak_100c);
plot(x100c_freq, x4c_spec);
plot(x100c_freq, x100c_spec);
legend('4C peaks', '100C peaks');

%xcorr_100c
%xcorr_50c
%xcorr_4c
% for i = 1:length(hp_spec)
%     for j = 1:length(locs1)
%         if hp_freq(i) == locs1(j)
%             peak_hp(i,1) = 1;
%         end
%     end
%     for k = 1:length(locs2)
%         if hp_freq(i) == locs2(k)
%             peak_100c(i,1) = 1;
%         end
%     end
%     for L = 1:length(locs1)
%         if hp_freq(i) == locs1(j)
%             peak_hp(i,1) = 1;
%         end
%     end
%     for m = 1:length(locs1)
%         if hp_freq(i) == locs1(j)
%             peak_hp(i,1) = 1;
%         end
%     end    
%end

            
pks1 = ones(length(pks1),1);
pks2 = ones(length(pks2),1);
pks3 = ones(length(pks3),1);
pks4 = ones(length(pks4),1);

% figure(1)
% plot(hp_freq,hp_spec)
% hold on 
% scatter(locs1,pks1)
% 
% figure(2)
% plot(x100c_freq,x100c_spec)
% hold on 
% scatter(locs2,pks2)
% 
% figure(3)
% plot(x50c_freq,x50c_spec)
% hold on 
% scatter(locs3,pks3)
% 
% figure(4)
% plot(x4c_freq,x4c_spec)
% hold on 
% scatter(locs4,pks4)
% 
% figure(5)
% scatter(locs1,pks1)
% hold on 
% scatter(locs2,pks2)
% hold on 
% scatter(locs3,pks3)
% hold on 
% scatter(locs4,pks4)

% figure(1)
% scatter(hp_freq,peak_hp)
% hold on
% scatter(hp_freq,peak_hp2)
% hold on
% plot(hp_freq,hp_spec)
% legend('jasons','peters','actual')

% figure(1)
% scatter(hp_freq,peak_hp)
% hold on
% scatter(hp_freq,peak_100c)
% hold on
% scatter(hp_freq,peak_50c)
% hold on
% scatter(hp_freq,peak_4c)
% legend('hp','100c','50c','4c')
