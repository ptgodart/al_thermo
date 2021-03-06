load('analysis.mat')
peak_hp = zeros(length(hp_spec),1);
peak_100c = zeros(length(hp_spec),1);
peak_50c = zeros(length(hp_spec),1);
peak_4c = zeros(length(hp_spec),1);

[pks1,locs1] = findpeaks(smooth(hp_spec,21),hp_freq);
[pks2,locs2] = findpeaks(smooth(x100c_spec,21),x100c_freq);
[pks3,locs3] = findpeaks(smooth(x50c_spec,21),x50c_freq);
[pks4,locs4] = findpeaks(smooth(x4c_spec,21),x4c_freq);
locs1 = [405, 501, 569, 615, 740, 885, 881, 964, 1570, 1454, 1165, 1074, 1400, 1631, 1975, 2104, 1811, 3097, 3303];
locs2 = [400, 478, 630, 731, 885, 950, 1016, 1068, 1159, 1400, 1518, 1635, 2092, 3103, 3408]; 
locs3 = [478, 526, 627, 732, 879, 975, 1034, 1072, 1159, 1398, 1531,1633, 2096, 3111, 3440];   
locs4 = [598, 858, 939, 1412, 1510, 1631, 2121, 3462];

hp_100c_score = 0;
hp_50c_score = 0;
hp_4c_score = 0;
score_4c_100c = 0;
score_4c_50c = 0;
peak_bound = 20;

for i = 1:length(locs1)
    for j = 1:length(locs2)
        if locs1(i) > (locs2(j)-peak_bound) && locs1(i) < (locs2(j)+peak_bound)
            hp_100c_score = hp_100c_score + 1 ;
        end
    end
    for k = 1:length(locs3)
        if locs1(i) > (locs3(k)-peak_bound) && locs1(i) < (locs3(k)+peak_bound)
            hp_50c_score = hp_50c_score + 1 ;
        end
    end
    for L = 1:length(locs4)
        if locs1(i) > (locs4(L)-peak_bound) && locs1(i) < (locs4(L)+peak_bound)
            hp_4c_score = hp_4c_score + 1 ;
        end
    end
end
hp_100c_score = hp_100c_score/length(locs2)
hp_50c_score = hp_50c_score/length(locs3)
hp_4c_score = hp_4c_score/length(locs4)


for i = 1:length(locs4)
    for m = 1:length(locs2)
        if locs4(i) > (locs2(m)-peak_bound) && locs4(i) < (locs2(m)+peak_bound)
            score_4c_100c = score_4c_100c + 1 ;
        end
    end
    for n = 1:length(locs3)
        if locs4(i) > (locs3(n)-peak_bound) && locs4(i) < (locs3(n)+peak_bound)
            score_4c_50c = score_4c_50c + 1 ;
        end
    end
end
score_4c_100c = score_4c_100c/length(locs2)
score_4c_50c = score_4c_100c/length(locs2)
% peak_hp(round((locs1'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_100c(round((locs2'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_50c(round((locs3'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;
% peak_4c(round((locs4'-hp_freq(1))/(hp_freq(2)-hp_freq(1))),:) = 1;

% xcorr_hp_100c = xcorr(peak_hp,peak_100c);
% figure(1)
% plot(xcorr_hp_100c)

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

            
% pks1 = ones(length(pks1),1);
% pks2 = ones(length(pks2),1);
% pks3 = ones(length(pks3),1);
% pks4 = ones(length(pks4),1);

figure(1)
plot(hp_freq,hp_spec)
hold on 
scatter(locs1,ones(length(locs1),1)*2)

figure(2)
plot(x100c_freq,x100c_spec)
hold on 
scatter(locs2,ones(length(locs2),1)*2)

figure(3)
plot(x50c_freq,x50c_spec)
hold on 
scatter(locs3,ones(length(locs3),1)*2)

figure(4)
plot(x4c_freq,x4c_spec)
hold on 
scatter(locs4,ones(length(locs4),1)*2)

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

<<<<<<< HEAD
% figure(1)
=======
% figure(5)
>>>>>>> 1df0478edf6b8bdae48868d8965e325100048a41
% scatter(hp_freq,peak_hp)
% hold on
% scatter(hp_freq,peak_100c)
% hold on
% scatter(hp_freq,peak_50c)
% hold on
% scatter(hp_freq,peak_4c)
% legend('hp','100c','50c','4c')
