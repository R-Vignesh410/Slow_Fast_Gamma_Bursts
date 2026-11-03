clc;clear
f=figure;
f.WindowState="Maximized";
plotHandles=getPlotHandles(1,4,[0.08 0.55 0.85 0.4],0.08,0.08,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
load(fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat'));
%badtrials file is loaded
load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
% selected electrode file is loaded-RF data file
load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
% LFP_data file is loaded
LFP_data_file=dir(fullfile(parent_file_path,"Decimated_8_LFP_data_alpa_H"));
load('gamma_duration_alpaH_MP.mat');
% natrisort the LFP data-remain in form of struct
LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'})); %remove . and ..
%I want all the files in the folder to be sorted by their name in order 1,2,3... (not 1,10,11,..)
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
%j=1;% counter for electrode position
electrode_num=length(selected_elec_LFP);
%%%%%%%%%%%%%%%%%%%%%%%%%%% LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 32- possible elec
selected_elec_pos=37;
i=selected_elec_LFP(selected_elec_pos);
load(fullfile('Decimated_8_LFP_Data_alpa_H',LFP_data_file{i}));%analogDataDecimatedDecimated
load('timeVals.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
current_SF=2; % 1 CPD
current_ORI=3; % 45 deg
trial=1; % first trial of SF- 1 cpd and 45 deg
% SF- 1 cpd and 157.5 deg- trial num=6, 13, 17
%%%%%%%%%%%%%% BAND-PASS FILTERING %%%%%%%%%%%%%%%%%%%%%%
trial_temp=[parameterCombinations{:,:,:,current_SF,current_ORI}];
trial_temp=setdiff(trial_temp,badTrials);
LFP_signal=analogDataDecimated(trial_temp,:);
% Plot the original LFP data signal
LFP_signal_one_trial=LFP_signal((trial),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresholdFraction=0.5;
num_iterations=120;
dict_Size=2500000;
jj=selected_elec_pos;
displayFlag=1;
% current elec- 41 (pos-37)
data_temp=LFP_signal_one_trial;
Fs=250;
params.Fs=Fs;
params.fpass=[0 150];
params.tapers=[1 1];
params.trialave=1;
params.pad=-1;
params.err=[2 0.05];
movingWin=[0.2 0.02];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S_temp,t_temp,f_temp]=mtspecgramc(data_temp',movingWin,params);
t_temp=t_temp+timeVals(1)-(1/Fs); % Center the times with respect to the stimulus onset time
subplot(plotHandles(1,2));
pcolor(t_temp,f_temp,log10(S_temp)'); hold on;
shading interp;
colormap('jet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,3));
[~,~]=getBurstLengthCGT(data_temp,timeVals,[],1,stimulusPeriodS,baselinePeriodS,[0 100],0.0125,2.5);
subplot(plotHandles(1,4));
[~,~]=getBurstLengthCGT(data_temp,timeVals,[],1,stimulusPeriodS,baselinePeriodS,[0 100],0.025,2.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gabor_SF=gaborInfo_accumulator_alpaH{current_SF,current_ORI,jj};
gabor_SF_one_trial=gabor_SF(trial,:,:);
header_SF=header_accumulator_alpaH{current_SF,current_ORI,jj};
header_SF_one_trial=header_SF(trial,:,:);
load('MP_non_filt_single_trial_results.mat')
%%%%%%%%%%%%%%%%%% Slow gamma %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_freq=[20 65];
num_iterations=600;
diffPower=getChangeInPower_single(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,gamma_freq);
thresholdFactor=sqrt(thresholdFraction*diffPower);
[~,~,~,~,~,~]= getBurstLengthMP_plot(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,gamma_freq,num_iterations,0.9,dict_Size,gabor_temp_one_trial(:,1:200,:),header_temp_one_trial,plotHandles,1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% plot the burst length at time burst center and y-axis frequency_measured
load("Detected_burst_data_MP_sample_elec_fig_S3.mat","length_measured_sg","freq_measured_sg","time_measured_sg","length_measured_fg","freq_measured_fg","time_measured_fg");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,1));
plot_bursts_on_TF_plot(length_measured_sg,length_measured_fg,...
    time_measured_sg,time_measured_fg,freq_measured_sg,freq_measured_fg,trial,...
    timeVals,slow_gamma_freq,fast_gamma_freq);
clim([-4 4]);
xlabel('Time(s)');
ylabel('Frequency(Hz)');
c=colorbar;
c.Position = c.Position + [0.03 0 0.003 0];
c.FontSize = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,2))
plot_bursts_on_TF_plot(length_measured_sg,length_measured_fg,...
    time_measured_sg,time_measured_fg,freq_measured_sg,freq_measured_fg,trial,...
    timeVals,slow_gamma_freq,fast_gamma_freq);
clim([-2 2]);
xlabel('Time(s)');
ylabel('Frequency(Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,3));
plot_bursts_on_TF_plot(length_measured_sg,length_measured_fg,...
    time_measured_sg,time_measured_fg,freq_measured_sg,freq_measured_fg,trial,...
    timeVals,slow_gamma_freq,fast_gamma_freq);
clim([3 9]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(plotHandles(1,4));
plot_bursts_on_TF_plot(length_measured_sg,length_measured_fg,...
    time_measured_sg,time_measured_fg,freq_measured_sg,freq_measured_fg,trial,...
    timeVals,slow_gamma_freq,fast_gamma_freq);
clim([3 9]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_axis_ticks_fontsize(plotHandles,22,18,1);
labels = {'A','B','C','D'};
x_positions = [0.03,0.27,0.505,0.74];
y_positions = [0.87];  
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%