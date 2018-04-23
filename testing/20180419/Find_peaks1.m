load('analysis.mat')
peak_hp = zeros(length(hp_spec),1);
peak_100c = zeros(length(hp_spec),1);
peak_50c = zeros(length(hp_spec),1);
peak_4c = zeros(length(hp_spec),1);

[pks1,locs1] = findpeaks(smooth(hp_spec,21),hp_freq);
[pks2,locs2] = findpeaks(smooth(x100c_spec,21),x100c_freq);
[pks3,locs3] = findpeaks(smooth(x50c_spec,21),x50c_freq);
[pks4,locs4] = findpeaks(smooth(x4c_spec,21),x4c_freq);

hp_100c_score = 0;
hp_50c_score = 0;
hp_4c_score = 0;
peak_bound = 15;

for i = 1:length(locs1)
    for j = 1:length(locs2)
        if locs1(i) > (locs2(j)-peak_bound) && locs1(i) < (locs2(j)+peak_bound)
            hp_100c_score = hp_100c_score + 1 ;
        end
    end
end
hp_100c_score = hp_100c_score/length(locs2)

for i = 1:length(locs1)
    for j = 1:length(locs3)
        if locs1(i) > (locs3(j)-peak_bound) && locs1(i) < (locs3(j)+peak_bound)
            hp_50c_score = hp_50c_score + 1 ;
        end
    end
end
hp_50c_score = hp_50c_score/length(locs3)

for i = 1:length(locs1)
    for j = 1:length(locs4)
        if locs1(i) > (locs4(j)-peak_bound) && locs1(i) < (locs4(j)+peak_bound)
            hp_4c_score = hp_4c_score + 1 ;
        end
    end
end
hp_4c_score = hp_4c_score/length(locs4)


% peak_hp(round((locs1'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_100c(round((locs2'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_50c(round((locs3'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_4c(round((locs4'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;

% xcorr_hp_100c = xcorr(peak_hp,peak_100c);
% figure(1)
% plot(xcorr_hp_100c)

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

            
% pks1 = ones(length(pks1),1);
% pks2 = ones(length(pks2),1);
% pks3 = ones(length(pks3),1);
% pks4 = ones(length(pks4),1);

figure(1)
plot(hp_freq,hp_spec)
hold on 
scatter(locs1,pks1)

figure(2)
plot(x100c_freq,x100c_spec)
hold on 
scatter(locs2,pks2)

figure(3)
plot(x50c_freq,x50c_spec)
hold on 
scatter(locs3,pks3)

figure(4)
plot(x4c_freq,x4c_spec)
hold on 
scatter(locs4,pks4)

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

% figure(5)
% scatter(hp_freq,peak_hp)
% hold on
% scatter(hp_freq,peak_100c)
% hold on
% scatter(hp_freq,peak_50c)
% hold on
% scatter(hp_freq,peak_4c)
% legend('hp','100c','50c','4c')
