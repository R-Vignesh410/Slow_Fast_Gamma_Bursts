f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(2,1,[0.48 0.52 0.44 0.42],0.04,0.05,0);
plotHandles_2=getPlotHandles(1,1,[0.48 0.08 0.44 0.36],0.04,0.08,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
load(fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat'));
% badtrials file is loaded
load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
% selected electrode file is loaded-RF data file
load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
% LFP_data file is loaded
LFP_data_file=dir(fullfile(parent_file_path,"Decimated_8_LFP_data_alpa_H"));
load('gamma_duration_alpaH_MP.mat');
% natrisort the LFP data-remain in form of struct
LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'})); %remove . and ..
% I want all the files in the folder to be sorted by their name in order 1,2,3... (not 1,10,11,..)
LFP_data_file = natsortfiles({LFP_data_file.name});
stimulusPeriodS=[0.25 0.75];
baselinePeriodS=[-0.5 0];
%%%%%%%% Used for calculating the change in power %%%%%%%%%%%%%%
fast_gamma_freq=[36 65];
%%%%%%% for alpaH - slow gamma range is taken from 25 to 40 Hz%%%%%%%%%
slow_gamma_freq=[20 32];
selected_elec=highRMSElectrodes;
selected_elec_LFP=selected_elec(1:77);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% j=1;% counter for electrode position
electrode_num=length(selected_elec_LFP);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selected_elec_pos=1;
i=selected_elec_LFP(selected_elec_pos);
load(fullfile('Decimated_8_LFP_Data_alpa_H',LFP_data_file{i})); % analogDataDecimatedDecimated
load('timeVals.mat');
LFP_non_decimated_path = fullfile( ...
    'C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\data\alpaH\Microelectrode\210817\GRF_002\segmentedData\LFP', ...
    sprintf('elec%d.mat', i));
load(LFP_non_decimated_path);
timeVals_non_decimated_path="C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab\data\alpaH\Microelectrode\210817\GRF_002\segmentedData\LFP\lfpInfo.mat";
Fs=250;
Fs_non_deci=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw dotted orange lines at 20 and 32 Hz and dotted blue lines at 36 and 65 Hz
color_orange=[0.9, 0.4, 0.0];
color_blue=[0 0 1];
current_SF=2; % 1 CPD
current_ORI=4; % 67.5 deg
trial=21;
%%%%%%%%%%%%%% BAND-PASS FILTERING %%%%%%%%%%%%%%%%%%%%%%
trial_temp=[parameterCombinations{:,:,:,current_SF,current_ORI}];
trial_temp=setdiff(trial_temp,badTrials);
LFP_signal=analogDataDecimated(trial_temp,:);
LFP_signal_non_deci=analogData(trial_temp,:);
% Plot the original LFP data signal
LFP_signal_one_trial=LFP_signal((trial),:);
LFP_signal_one_trial_non_deci=LFP_signal_non_deci((trial),:);
% filter between 20 Hz and 32 Hz
slow_gamma_freq_burst_range=[20 32];
[b,a]=butter(4,slow_gamma_freq_burst_range/(Fs_non_deci/2),'bandpass');
LFP_signal_one_trial_sg=filtfilt(b,a,LFP_signal_one_trial_non_deci);
% filter between 36 Hz and 65 Hz
fast_gamma_freq_burst_range=[36 65];
% fast_gamma_freq_burst_range=[40 65];
[b,a]=butter(4,fast_gamma_freq_burst_range/(Fs_non_deci/2),'bandpass');
LFP_signal_one_trial_fg=filtfilt(b,a,LFP_signal_one_trial_non_deci);
% filter between 3 and 100 Hz
wide_band_freq=[3 100];
[b,a]=butter(4,wide_band_freq/(Fs_non_deci/2),'bandpass');
LFP_signal_one_trial_wb=filtfilt(b,a,LFP_signal_one_trial_non_deci);
% plot the signal
load(timeVals_non_decimated_path);
%%%%%%%%%%%%%%%%%%%%%%% SLOW GAMMA %%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(2,1));
% plot(timeVals,LFP_signal_one_trial_sg,'-b','LineWidth',2.5);
plot(timeVals,LFP_signal_one_trial_sg,'-b','LineWidth',1.25);
hold on;
plot(timeVals,LFP_signal_one_trial_wb,'-k','LineWidth',1.25);
ylim([-150 150]);
%%%%%%%%%%%%%%%%%%%%% FAST GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(1,1));
% plot(timeVals,LFP_signal_one_trial_fg,'-','Color',color_orange,'LineWidth',2.5);
plot(timeVals,LFP_signal_one_trial_fg,'-','Color',color_orange,'LineWidth',1.25);
hold on;
plot(timeVals,LFP_signal_one_trial_wb,'-k','LineWidth',1.25);
ylabel("Amplitude (\muV)");
ylim([-150 150]);
%%%%%%%%%%%%%%%%%%%% LFP signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("timeVals.mat");
thresholdFraction=0.25;
% thresholdFraction=0.5;
num_iterations=120;
dict_Size=2500000;
jj=selected_elec_pos;
displayFlag=1;
data_temp=LFP_signal_one_trial;
Fs=250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gabor_SF=gaborInfo_accumulator_alpaH{current_SF,current_ORI,jj};
gabor_SF_one_trial=gabor_SF(trial,:,:);
header_SF=header_accumulator_alpaH{current_SF,current_ORI,jj};
header_SF_one_trial=header_SF(trial,:,:);
%%%%%%%%%%%%%%%%%% Slow gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_freq=[20 65];
diffPower=getChangeInPower_single(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,gamma_freq);
thresholdFactor=sqrt(thresholdFraction*diffPower);
idx_1=1;idx_2=1;
[length_measured_broad,freq_measured_broad,time_measured_broad,~,~,~]= getBurstLengthMP_plot(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_Size,gabor_SF_one_trial,header_SF_one_trial,plotHandles_4,idx_1,idx_2);
subplot(plotHandles_4(1,1));
hold on;
% set y-axis to [0 100] and x axis to the timeVals
axis([timeVals(1) timeVals(end) 0 100]);
% plot the  horizontal line (Across time Vals)- frequency ranges using dashed white lines(sloww gamma) and dashed black lines(fast gamma)
line([timeVals(1) timeVals(end)],[slow_gamma_freq(1) slow_gamma_freq(1)],'Color','k','LineStyle','--','LineWidth',2);
line([timeVals(1) timeVals(end)],[slow_gamma_freq(2) slow_gamma_freq(2)],'Color','k','LineStyle','--','LineWidth',2);
% plot the  horizontal line (Across time Vals)- frequency ranges using dashed white lines(sloww gamma) and dashed black lines(fast gamma)
line([timeVals(1) timeVals(end)],[fast_gamma_freq(1) fast_gamma_freq(1)],'Color','k','LineStyle','--','LineWidth',2);
line([timeVals(1) timeVals(end)],[fast_gamma_freq(2) fast_gamma_freq(2)],'Color','k','LineStyle','--','LineWidth',2);
% plot the burst length at time burst center and y-axis frequency_measured
load('Detected_bursts_MP_example_elec_fig_1.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_3(1,1));
hold on;
plot_burst_lengths(length_measured_fg,time_measured_fg,trial,'m');
subplot(plotHandles_3(2,1));
hold on;
plot_burst_lengths(length_measured_sg,time_measured_sg,trial,'g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles_4(1,1));
hold on;
for ii = 1:length(length_measured_sg{trial})
    if length_measured_sg{trial}(ii)>0.8
        continue
    end
    start_time = time_measured_sg{trial}(ii) - (length_measured_sg{trial}(ii)/2);
    end_time = time_measured_sg{trial}(ii) + (length_measured_sg{trial}(ii)/2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot([start_time,end_time],[freq_measured_sg{trial}(ii),freq_measured_sg{trial}(ii)],'Color','k','LineWidth',2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot(time_measured_sg{trial}(ii),freq_measured_sg{trial}(ii),'ko','MarkerSize',7,'MarkerFaceColor','g');
end
for ii=1:length(length_measured_fg{trial})
    if length_measured_fg{trial}(ii)>0.8
        continue
    end
    start_time = time_measured_fg{trial}(ii) - length_measured_fg{trial}(ii)/2;
    end_time = time_measured_fg{trial}(ii) + length_measured_fg{trial}(ii)/2;
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot([start_time,end_time],[freq_measured_fg{trial}(ii),freq_measured_fg{trial}(ii)],'Color','k','LineWidth',2);
    % plot the burst length at time burst center- same as before and y-axis frequency_measured
    plot(time_measured_fg{trial}(ii),freq_measured_fg{trial}(ii),'ko','MarkerSize',7,'MarkerFaceColor','m');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % draw a white dotted line at 0.25 and 0.75
line([0.25 0.25], [0 100], 'Color','m', 'LineStyle', '--','LineWidth',2);
line([0.75 0.75], [0 100], 'Color', 'm', 'LineStyle', '--','LineWidth',2);
line([-0.5 -0.5], [-100 100], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
% line([0 0], [-100 100], 'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
% line([0.8 0.8],[-100 100],'Color',[0.23,0.51,0], 'LineStyle', '--','LineWidth',2);
line([0 0], [-100 100], 'Color','w', 'LineStyle', '--','LineWidth',2);
line([0.8 0.8],[-100 100],'Color','w', 'LineStyle', '--','LineWidth',2);
c=colorbar;
c.Position = c.Position + [ 0.04 0 0.005 0];
c.Ticks  = [-4 -2 0 2 4];
c.FontSize = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis([-0.1 1 0 80]);
clim([-4 4]);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_axis_ticks_fontsize(plotHandles_3,22,18,1);
set_axis_ticks_fontsize(plotHandles_3,22,18,2);
set_axis_ticks_fontsize(plotHandles_4,22,18,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=annotation('textbox',...
[0.93 0.575 0.06 0.06],...
'String',{'Slow \gamma'},...
'Rotation',0,...
'FontWeight','bold',...
'FontSize',22,...
'Color','b',...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
t1.FitBoxToText = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2=annotation('textbox',...
[0.93 0.81 0.06 0.06],...
'String',{'Fast \gamma'},...
'Rotation',0,...
'FontWeight','bold',...
'FontSize',22,...
'Color',color_orange,...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
t2.FitBoxToText = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t3=annotation('textbox',...
[0.985 0.085 0.06 0.06],...
'String',{'Raw Power (log-scale)'},...
'Rotation',90,...
'FontSize',22,...
'Color','k',...
'FontName','Helvetica',...
'EdgeColor',[1 1 1]);
t3.FitBoxToText = 'on';
labels = {'A','B'};
x_positions = [0.425];
y_positions = [0.85,0.365];  
k = 1;
for j = 1:length(y_positions)
    for i = 1:length(x_positions)
        annotation('textbox', ...
            [x_positions(i), y_positions(j), 0.1, 0.1], ...
            'String', labels{k}, ...
            'FontSize', 28, ...
            'FontWeight', 'Bold', ...
            'EdgeColor', 'none', ...
            'FontName', 'Helvetica');
        k = k + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%