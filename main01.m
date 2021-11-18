clear all; clear clc; close all;
global R; %R peaks vector
is_noised = 0; % 0 - clean, 1 - noised signal

[filename,folder] = uigetfile('*');
path = fullfile(folder,filename);
signal = importdata(path);

f_s = 1000; %Hz
set_fs(f_s); %Hz
T=1/f_s; %sec
time_vec = 0:T:(length(signal)-1)*T;
set_time_vec(time_vec);

if is_noised == 1
    set_sig_filtered(total_filter(signal, f_s, 2, 50,100, 40));
    normalized_signal=normalized(get_sig_filtered());  %normalizing the filter
else
    set_sig_filtered(signal);
    normalized_signal=normalized(signal);
end
%pre-QRS detecion
linear_filtered_butterworth=linearFILTER(normalized_signal,10,25,f_s); %linear filtering
amplified_ECG_signal = amplified(linear_filtered_butterworth.*1000,3); %non-linear transformation
set_amplified_ECG_signal(amplified_ECG_signal); %non-linear transformation

figure;
uicontrol('style','push',...
    'units','pix',...
    'position',[400 5 100 20],...
    'fontsize',8,...
    'string','Set Threshold',...
    'callback',{@pb_call});

plot(time_vec,amplified_ECG_signal); % creates a figure and plot something on it
title('amplified ECG signal- non linear transformation ');
xlabel('time [sec]');
ylabel('arbitrary units');

%creating input box
msgbox('in order to set a threshold, get close to the signal untill you can determine a value that is greater than most T waves. push the "Set Threshold" button when you ready', 'Instructions');

function pb_call(src, event)
f_s=get_fs();
sig_filtered = get_sig_filtered();
time_vec = get_time_vec();
amplified_ECG_signal = get_amplified_ECG_signal();
prompt = {'Enter threshold value'};
dlgtitle = 'threshold value';
definput = {'100'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,[1 40],definput,opts);
user_threshold=str2num(answer{1});
figure;
plot(time_vec,sig_filtered);
hold on
%plot(R_peaks(:,4),R_peaks(:,3),'r*');
%}
R = findRpeaks(sig_filtered, amplified_ECG_signal,user_threshold, 90, 20);
plot(R./1000,sig_filtered(R),'r*');
title('the original signal with the R waves marked');
xlabel('[sec]');
ylabel('[mV]');

%ploting HRV
peakInterval = diff(R);
locs_Rwave_time=(R-1)/f_s;
TimeInterval=peakInterval/f_s; %calculating time intervals between peaks
BPM_ECG=1./(TimeInterval)*60; %calculatin BPM over time
BPM_ECG_time=locs_Rwave_time(2:length(locs_Rwave_time));
figure
plot(BPM_ECG_time, BPM_ECG);  %plotin the ECG pulse
xlabel('time [s]');
ylabel('Beats');
title('ECG - HRV');
set_R(R);
end


function set_R(val)
global R
R = val;
end
function r = get_R
global R
r = R;
end

function set_fs(val)
global fs
fs = val;
end
function r = get_fs
global fs
r = fs;
end

function set_time_vec(val)
global timevec
timevec = val;
end
function r = get_time_vec
global timevec
r = timevec;
end

function set_sig_filtered(val)
global sig_filtered
sig_filtered = val;
end
function r = get_sig_filtered
global sig_filtered
r = sig_filtered;
end

function set_amplified_ECG_signal(val)
global amplified_ECG_signal
amplified_ECG_signal = val;
end
function r = get_amplified_ECG_signal
global amplified_ECG_signal
r = amplified_ECG_signal;
end








% R_peaks = run_methods2(signal, is_noised);

%-------------------------------functions-------------------------------

function R = run_methods(signal, is_noised)
f_s=1000; %Hz
T=1/f_s; %sec
time_vec = 0:T:(length(signal)-1)*T;

if is_noised == 1
    sig_filtered=total_filter(signal, f_s, 2, 50,100, 40);
    normalized_signal=normalized(sig_filtered);  %normalizing the filter
else
    sig_filtered = signal;
    normalized_signal=normalized(signal);
end
%pre-QRS detecion
linear_filtered_butterworth=linearFILTER(normalized_signal,10,25,f_s); %linear filtering
amplified_ECG_signal=amplified(linear_filtered_butterworth.*1000,3); %non-linear transformation


figure
plot(time_vec,amplified_ECG_signal);
title('amplified ECG signal- non linear transformation ');
xlabel('time [sec]');
ylabel('arbitrary units');
%creating input box
f=msgbox('in order to set a threshold- get close to the signal- untill you can determine a value that is greater than most T waves', 'pre- threshold setting');
%creating input box
f=1;
prompt = {'Enter threshold value in range 60-250 acording to the previos figure'};
dlgtitle = 'threshold value';
definput = {'100'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,[1 40],definput,opts);
user_threshold=str2num(answer{1});


figure;
plot(time_vec,sig_filtered);
hold on
%plot(R_peaks(:,4),R_peaks(:,3),'r*');
%}
R = findRpeaks(sig_filtered, amplified_ECG_signal,user_threshold, 90, 20);
plot(R./1000,sig_filtered(R),'r*');
title('the original signal with the R waves marked');
xlabel('[sec]');
ylabel('[mV]');

%ploting HRV
peakInterval = diff(R);
locs_Rwave_time=(R-1)/f_s;
TimeInterval=peakInterval/f_s; %calculating time intervals between peaks
BPM_ECG=1./(TimeInterval)*60; %calculatin BPM over time
BPM_ECG_time=locs_Rwave_time(2:length(locs_Rwave_time));
figure
plot(BPM_ECG_time, BPM_ECG);  %plotin the ECG pulse
xlabel('time [s]');
ylabel('Beats');
title('ECG - HRV');
end

function normalized_signal=normalized(ECG_signal)
%normalizing the input signal
max_value=max(ECG_signal);
min_value=min(ECG_signal);
normalized_signal=-1+((ECG_signal-min_value)./(max_value-min_value)).*2;
end

function linear_filtered_signal=linearFILTER(ECG_signal,fl,fh,fs)
%butterworth high and low pass filters at f_low Hz and f_high Hz = BPF
%annuating the P and T waves - in our case fl=10 Hz, fh=25 Hz
wn_H=(fh/(fs/2))-0.002; % f/(fs/2)
wn_L=(fl/(fs/2))+0.002; % f/(fs/2)
[b_high,a_high] = butter(5,wn_H,'high');
[b_low,a_low] = butter(5,wn_L,'low');
filtered_high_ECG=filter(b_high,a_high,ECG_signal);
filtered_low_ECG=filter(b_low,a_low,filtered_high_ECG);
linear_filtered_signal=filtered_low_ECG;
end

function amplified_ECG_signal=amplified(ECG_signal,alpha)
%non linear transformation - amplifiying the QRS complex
amplified_ECG_signal=ECG_signal.^alpha;
end

function ambient_noise_remove=specipic_Hz_filter(ECG_signal,f1,f2,fs)
%using butterworth filter to remove 50 Hz and 100 Hz
wn1=[f1/(fs/2)-0.002,f1/(fs/2)+0.002];% f/(fs/2)
[b1,a1]=butter(2,wn1,'stop'); %butterworth filter
wn2=[f2/(fs/2)-0.002,f2/(fs/2)+0.002];
[b2,a2]=butter(2,wn1,'stop');
filtered_f1Hz=filter(b1,a1,ECG_signal);
ambient_noise_remove=filter(b2,a2,filtered_f1Hz);
end

function signal_removed_baseline_wander=remove_baseline(ECG_signal,fc,fs)
%removing the baseline wander - in our case cutoff frequency 2 Hz
HPF=fir1(fs,fc/(fs/2),'high');% apply FIR HPF.
signal_removed_baseline_wander=filter(HPF,1,ECG_signal);
end

function Rpeaks = findRpeaks(ECG_orig, ECG_filt, threshold, win_left, win_right)
%finding the R peaks at the original signal
%threshold=100;
ECG_after_threshold = ECG_filt>threshold;
R_peaks = zeros(1000,4);
n = 1;
first1index=find(ECG_after_threshold,1);
len = length(ECG_filt);
while first1index<len
    R_peaks(n,1) = first1index;
    last1index = find(ECG_after_threshold(first1index:end)==0,1)+first1index-2;
    if isempty(last1index)
        last1index=len;
    end
    R_peaks(n,2) = last1index;
    [R_peaks(n,3), R_peaks(n,4)] = max(ECG_orig(max(1,first1index-win_left):min(last1index+win_right,len)));
    R_peaks(n,4) = R_peaks(n,4) + first1index-win_left -1;
    first1index = find(ECG_after_threshold(last1index+1:end),1) + last1index;
    n = n+1;
end
Rpeaks = R_peaks(1:n-1,4);
end


function signal_LPF=LPF_first(ECG_signal, f_s, fc)
wn=(fc/(f_s/2));
[b,a]=butter(6,wn);
signal_LPF=filter(b,a,ECG_signal);
end


function filtered_signal=total_filter(ECG_signal, f_s, basline_wander_freq, raash_reshet1,raash_reshet2, fc)
filtered15=LPF_first(ECG_signal, f_s, fc);
filtered2=specipic_Hz_filter(filtered15,raash_reshet1,raash_reshet2,f_s); %%using butterworth filter to remove 50 Hz and 100 Hz
filtered3=remove_baseline(filtered2,basline_wander_freq,f_s);% %removing the baseline wander - in our case cutoff frequency 2 Hz
filtered_signal=filtered3;
end