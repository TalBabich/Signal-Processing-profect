clear all; close all;
peaks = load('Peaks.mat');
R = cell(1,10);
for i=1:10
    R{i} = peaks.Peaks{i,1}.PeakDetection.Rwaves.';
end

% the rest window in every signal
rest_time = [20,550;
    880,1060;
    1200,1700;
    1,330;
    840,1090;
    10,309;
    10,400;
    10,432;
    70,1100;
    1473,1815];

% the ex window in every signal
ex_time = [900,1400;
    480,690;
    450,920;
    600,860;
    330,600;
    426,730;
    1050,1640;
    830,1140;
    1550,2283;
    680,1350];

%definig time and frequency veriabels
f_s=1000; %Hz
T_input=1/f_s; %sec
fr=6; %Hz

% for i=1:10
%     figure(i);
%     RR = RR_calculation(R{i});
%     plot((1000*60)./RR);
% end

berger_rest = cell(2,10); %1 - berger values, 2 - time vector
berger_ex = cell(2,10);
pxx_rest = cell(2,10); %1 - pxx values, 2 - w
pxx_ex = cell(2,10);
ratio_rest = zeros(1,10); %LF/HF
ratio_ex = zeros(1,10);
for i=1:10
    R_rest = R{i}(rest_time(i,1):rest_time(i,2));
    R_ex = R{i}(ex_time(i,1):ex_time(i,2));
    [R_rest_fixed, NN_rest] = NN_calc(R_rest);
    [R_ex_fixed, NN_ex] = NN_calc(R_ex);
    %     figure(i);
    %     plot((1000*60)./NN_rest, 'r');
    %     hold on
    %     plot((1000*60)./NN_ex, 'b');
    [berger_rest{1,i}, berger_rest{2,i}] = berger(fr, R_rest_fixed, NN_rest);
    [berger_ex{1,i}, berger_ex{2,i}] = berger(fr, R_ex_fixed, NN_ex);
    [pxx_rest{1,i}, pxx_rest{2,i}] = pwelch(berger_rest{1,i});
    pxx_rest{2,i} = pxx_rest{2,i}./(2*pi);
    [pxx_ex{1,i}, pxx_ex{2,i}] = pwelch(berger_ex{1,i});
    pxx_ex{2,i} = pxx_ex{2,i}./(2*pi);
    LF = sum(pxx_rest{1,i}(pxx_rest{2,i}>0.04 & pxx_rest{2,i}<0.15));
    HF = sum(pxx_rest{1,i}(pxx_rest{2,i}>0.15 & pxx_rest{2,i}<0.4));
    ratio_rest(i) = LF/HF;
    LF = sum(pxx_ex{1,i}(pxx_ex{2,i}>0.04 & pxx_ex{2,i}<0.15));
    HF = sum(pxx_ex{1,i}(pxx_ex{2,i}>0.15 & pxx_ex{2,i}<0.4));
    ratio_ex(i) = LF/HF;
end

%signals 8-10 plotting
figure;
plot(berger_rest{2,8},60.*berger_rest{1,8}, 'r');
hold on
plot(berger_rest{2,9},60.*berger_rest{1,9}, 'b');
hold on
plot(berger_rest{2,10},60.*berger_rest{1,10}, 'g');
hold on
plot(berger_ex{2,8},60.*berger_ex{1,8}, 'r');
hold on
plot(berger_ex{2,9},60.*berger_ex{1,9}, 'b');
hold on
plot(berger_ex{2,10},60.*berger_ex{1,10}, 'g');
title('HR vs Time')
ylabel('HR');
xlabel('time [sec]');
legend('signal #8', 'signal #9', 'signal #10');

%signals 8-10 PSD plotting
for i=8:10
    figure;
    plot(pxx_rest{2,i}, 10.*log(pxx_rest{1,i}),'DisplayName','rest');
    hold on
    plot(pxx_ex{2,i}, 10.*log(pxx_ex{1,i}),'DisplayName','exercise');
    title(['PSD - signal #', num2str(i)]);
    xlabel('f [Hz]');
    ylabel('PSD');
    legend;
end

%errorbar plotting
figure;
x = categorical({'rest','exersice'});
bar(x, [mean(ratio_rest), mean(ratio_ex)], 'FaceColor', '#4DBEEE', 'EdgeColor','blue','LineWidth',2);
hold on
errorbar(x, [mean(ratio_rest), mean(ratio_ex)], [std(ratio_rest), std(ratio_ex)], 'o', 'LineWidth', 2, 'Color', '#D95319');
ylabel('Ratio LF/HF');

%Get the R peaks vector
%Return the NN vector and new R vector
function [R_fixed, NN] = NN_calc(R_vec)
NN = zeros(1, length(R_vec)-1);
for i=1:length(R_vec)-1
    NN(i)=R_vec(i+1)-R_vec(i);
end
avg = mean(NN(:));
NN(find(NN<0.8*avg | NN>1.2*avg)) = avg;
R_fixed = zeros(1, length(R_vec));
R_fixed(1) = R_vec(1);
for i=1:length(R_vec)-1
    R_fixed(i+1) = R_fixed(i) + NN(i);
end
end

%Get fr, R peaks vector and NN vector
%Return r vector after berger algorithm and a time vector (when every r value was measured).
function [r_vec, time_vec] = berger(fr, R_vec, NN_vec)
Tr = floor(1000/fr); %msec
r_vec = zeros(1, floor((R_vec(end)-R_vec(1))/Tr)-2);
time_vec = r_vec;
k = 1;
i = ceil(R_vec(1)/Tr)*Tr;
while i <= R_vec(end)-2*Tr
    Rindex = find(R_vec>i & R_vec<i+2*Tr, 1);
    if isempty(Rindex)
        r_vec(k) = (2*Tr)/NN_vec(find(R_vec>i,1)-1);
    else
        r_vec(k) = (i+2*Tr-R_vec(Rindex))/NN_vec(Rindex) + (R_vec(Rindex)-i)/NN_vec(Rindex-1);
    end
    i = i+Tr;
    time_vec(k) = i;
    k=k+1;
end
r_vec = r_vec.*(1000/(2*Tr));
end

% function RR_vec=RR_calculation(Rwaves_vector)
% %creating the RR intervals vector
% RR_vec=zeros(1,(length(Rwaves_vector)-1));
% for i=1:length(Rwaves_vector)-1
%     RR_vec(i)=Rwaves_vector(i+1)-Rwaves_vector(i);
% end
% end
%
% function [R_fixed, NN_fixed] = NN_movingAvg(Rpeaks, win)
% %creating Normal to Normal intervals
% RR = zeros(1, length(Rpeaks)-1);
% for i=1:length(Rpeaks)-1
%     RR(i)=Rpeaks(i+1)-Rpeaks(i);
% end
% NN_fixed = RR;
% for i=(win+1):(length(RR)-win)
%     avg = mean(RR(i-win:i+win));
%     if NN_fixed(i)<0.8*avg || NN_fixed(i)>1.2*avg
%         NN_fixed(i) = avg;
%     end
% end
% R_fixed = zeros(1, length(Rpeaks));
% R_fixed(1) = Rpeaks(1);
% for i=1:length(Rpeaks)-1
%     R_fixed(i+1) = R_fixed(i) + NN_fixed(i);
% end
% end
