function plot_overlap_sg_fg(plotHandles,Monkey_num)
    parent_file_path='C:\Users\rviiy\OneDrive - Indian Institute of Science\gamma_length_project_EEG_SRAYlab';
    displayFlag=0;
    stimulusPeriodS=[0.25 0.75];
    baselinePeriodS=[-0.5 0];
    thresholdFraction=0.5;
    num_iterations=120;
    dict_size=2500000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Monkey_num==1
        % Monkey- alpaH
        load('gamma_duration_alpaH_MP.mat')
        load((fullfile(parent_file_path,'alpaH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'alpaH_info','badTrials.mat'));
        load(fullfile(parent_file_path,'alpaH_info','alpaHMicroelectrodeRFData.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 32];
        fast_gamma_freq=[36 65];
        gabor_accumulator=gaborInfo_accumulator_alpaH;
        header_accumulator=header_accumulator_alpaH;
    elseif Monkey_num==2
        % Monkey- kesariH
        load('gamma_duration_kesariH_MP.mat')
        load((fullfile(parent_file_path,'kesariH_info','parameterCombinations.mat')))
        load(fullfile(parent_file_path,'kesariH_info','badTrials_kesari.mat'));
        load(fullfile(parent_file_path,'kesariH_info','kesariHMicroelectrodeRFData_Two.mat'));
        LFP_data_file=dir(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H'));
        LFP_data_file = LFP_data_file(~ismember({LFP_data_file.name},{'.','..'}));
        LFP_data_file = natsortfiles({LFP_data_file.name});
        slow_gamma_freq=[20 38];
        fast_gamma_freq=[42 65];
        gabor_accumulator=gaborInfo_accumulator_kesariH;
        header_accumulator=header_accumulator_kesariH;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('timeVals.mat')
    num_elec=length(current_electrode);
    counter=1;
    total_trial_counter=0;
    sg_only_trial_counter=0;
    fg_only_trial_counter=0;
    sg_fg_trial_counter=0;
    sg_fg_overlap_trial_counter=0;
    sg_bursts_occupancy_gatherer=cell(1,num_elec);
    fg_bursts_occupancy_gatherer=cell(1,num_elec);
    sg_fg_overlap_occupancy_gatherer=cell(1,num_elec);
    sg_fg_overlap_occupancy_gatherer_shuffle_trials=cell(1,num_elec);
    length_gatherer_sg=cell(1,num_elec);
    length_gatherer_fg=cell(1,num_elec);
    power_gatherer_sg=zeros(1,num_elec);
    power_gatherer_fg=zeros(1,num_elec);
    rng(2); % Set seed for reproducibility
    for i=current_electrode
        if Monkey_num==1
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_alpa_H',LFP_data_file{i}))
        end
        if Monkey_num==2
            load(fullfile(parent_file_path,'Decimated_8_LFP_data_kesari_H',LFP_data_file{i}))
        end
        % ORI- 157.5 deg and SF- 1cpd : for analysis
        ORI_num=8; SF_num=2;
        trial_temp=parameterCombinations{:,:,:,SF_num,ORI_num};
        trial_temp=setdiff(trial_temp,badTrials);
        data_temp=analogDataDecimated(trial_temp,:);
        gabor_temp=gabor_accumulator{SF_num,ORI_num,counter};
        header_temp=header_accumulator{SF_num,ORI_num,counter};
        %%%%%% Slow gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,slow_gamma_freq);
        power_gatherer_sg(counter)=diffPower;
        thresholdFactor=thresholdFraction*sqrt(diffPower);
        [length_temp_sg,freq_temp_sg,time_center_temp_sg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,slow_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        %%%%%% Fast gamma burst computation %%%%%%
        diffPower=getChangeInPower(data_temp,timeVals,stimulusPeriodS,baselinePeriodS,fast_gamma_freq);
        power_gatherer_fg(counter)=diffPower;
        thresholdFactor=thresholdFraction*sqrt(diffPower);
        [length_temp_fg,freq_temp_fg,time_center_temp_fg,~,~,~]= getBurstLengthMP(data_temp,timeVals,thresholdFactor,displayFlag,stimulusPeriodS,baselinePeriodS,fast_gamma_freq,num_iterations,0.9,dict_size,gabor_temp,header_temp);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        length_temp_sg_all_trials=[];
        length_temp_fg_all_trials=[];
        for ii=1:length(length_temp_sg)
            if ~isempty(length_temp_sg{ii})
                reject_idx_sg=find((length_temp_sg{ii}')>0.8);
                length_temp_sg{ii}(reject_idx_sg)=[];
                freq_temp_sg{ii}(reject_idx_sg)=[];
                time_center_temp_sg{ii}(reject_idx_sg)=[];
            end
            length_temp_sg_all_trials=[length_temp_sg_all_trials, length_temp_sg{ii}'];
        end
        for ii=1:length(length_temp_fg)
            if ~isempty(length_temp_fg{ii})
                reject_idx_fg=find((length_temp_fg{ii}')>0.8);
                length_temp_fg{ii}(reject_idx_fg)=[];
                freq_temp_fg{ii}(reject_idx_fg)=[];
                time_center_temp_fg{ii}(reject_idx_fg)=[];
            end
            length_temp_fg_all_trials=[length_temp_fg_all_trials, length_temp_fg{ii}'];
        end
        % arrange 1:length(length_temp_sg) randomly (shuffled)
        shuffled_indices=randperm(length(length_temp_sg));
        for ii=1:length(length_temp_sg)
            total_trial_counter=total_trial_counter+1;
            if (isempty(length_temp_sg{ii})||isempty(length_temp_fg{ii}))
                if ~isempty(length_temp_sg{ii})
                    sg_only_trial_counter=sg_only_trial_counter+1;
                    [sg_burts_union_length, ~, ~, ~, ~] = get_overlap_bursts(time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5),time_center_temp_sg{ii}'+((length_temp_sg{ii}')*0.5),[],[]);
                     sg_bursts_occupancy_gatherer{counter}=[sg_bursts_occupancy_gatherer{counter}; (sg_burts_union_length)];
                end
                if ~isempty(length_temp_fg{ii})
                    fg_only_trial_counter=fg_only_trial_counter+1;
                    [~,fg_burts_union_length, ~, ~, ~] = get_overlap_bursts([],[],time_center_temp_fg{ii}'-((length_temp_fg{ii}')*0.5),time_center_temp_fg{ii}'+((length_temp_fg{ii}')*0.5));
                    fg_bursts_occupancy_gatherer{counter}=[fg_bursts_occupancy_gatherer{counter}; (fg_burts_union_length)];
                end
                continue;
            end
            if (~isempty(length_temp_sg{ii}) && ~isempty(length_temp_fg{ii}))
                sg_fg_trial_counter=sg_fg_trial_counter+1;
                % Check for overlap between slow and fast gamma bursts
                start_point_temp_sg=time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5);
                end_point_temp_sg=time_center_temp_sg{ii}'+((length_temp_sg{ii}')*0.5);
                start_point_temp_fg=time_center_temp_fg{ii}'-((length_temp_fg{ii}')*0.5);
                end_point_temp_fg=time_center_temp_fg{ii}'+((length_temp_fg{ii}')*0.5);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [sg_bursts_union_length, fg_bursts_union_length, sg_fg_overlap_length, ...
                ~, ~] = get_overlap_bursts(start_point_temp_sg,end_point_temp_sg,...
                start_point_temp_fg, end_point_temp_fg);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sg_bursts_occupancy_gatherer{counter}=[sg_bursts_occupancy_gatherer{counter}; (sg_bursts_union_length)];
                fg_bursts_occupancy_gatherer{counter}=[fg_bursts_occupancy_gatherer{counter}; (fg_bursts_union_length)];
                norm_fact=1;
                sg_fg_overlap_occupancy_gatherer{counter}=[sg_fg_overlap_occupancy_gatherer{counter}; (sg_fg_overlap_length/norm_fact)];
                overlap_flag=(sg_fg_overlap_length>0);
                if overlap_flag
                    sg_fg_overlap_trial_counter=sg_fg_overlap_trial_counter+1;
                end
            end
            if (~isempty(length_temp_sg{ii}) && ~isempty(length_temp_fg{shuffled_indices(ii)}))
                % Check for overlap between slow and fast gamma bursts
                start_point_temp_sg=time_center_temp_sg{ii}'-((length_temp_sg{ii}')*0.5);
                end_point_temp_sg=time_center_temp_sg{ii}'+((length_temp_sg{ii}')*0.5);
                start_point_temp_fg=time_center_temp_fg{shuffled_indices(ii)}'-((length_temp_fg{shuffled_indices(ii)}')*0.5);
                end_point_temp_fg=time_center_temp_fg{shuffled_indices(ii)}'+((length_temp_fg{shuffled_indices(ii)}')*0.5);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [sg_bursts_union_length, fg_bursts_union_length, sg_fg_overlap_length, ...
                 ~, ~] = get_overlap_bursts(start_point_temp_sg,...
                end_point_temp_sg, start_point_temp_fg, end_point_temp_fg);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                norm_fact=1;
                sg_fg_overlap_occupancy_gatherer_shuffle_trials{counter}=[sg_fg_overlap_occupancy_gatherer_shuffle_trials{counter}; (sg_fg_overlap_length/norm_fact)];
            end  
        end
        length_gatherer_sg{counter}=length_temp_sg_all_trials;
        length_gatherer_fg{counter}=length_temp_fg_all_trials;
        counter=counter+1;
    end
    % Display trial counts
    disp(['Total trials: ', num2str(total_trial_counter),'for Monkey_num: ', num2str(Monkey_num)]);
    disp(['Slow gamma only trials: ', num2str(sg_only_trial_counter),'for Monkey_num: ', num2str(Monkey_num)]);
    disp(['Fast gamma only trials: ', num2str(fg_only_trial_counter),'for Monkey_num: ', num2str(Monkey_num)]);
    disp(['Slow and fast gamma trials: ', num2str(sg_fg_trial_counter),'for Monkey_num: ', num2str(Monkey_num)]);
    disp(['Slow and fast gamma overlap trials: ', num2str(sg_fg_overlap_trial_counter),'for Monkey_num: ', num2str(Monkey_num)]);
    mean_overlap_sg_fg_gatherer=zeros(1,num_elec);
    mean_overlap_sg_fg_gatherer_shuffle_trials=zeros(1,num_elec);
    mean_burst_occupancy_sg_gatherer=zeros(1,num_elec);
    mean_burst_occupancy_fg_gatherer=zeros(1,num_elec);
    for i=1:num_elec
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mean_overlap_sg_fg_gatherer(i)=mean(sg_fg_overlap_occupancy_gatherer{i});
        mean_overlap_sg_fg_gatherer_shuffle_trials(i)=mean(sg_fg_overlap_occupancy_gatherer_shuffle_trials{i});
        mean_burst_occupancy_sg_gatherer(i)=mean(sg_bursts_occupancy_gatherer{i});
        mean_burst_occupancy_fg_gatherer(i)=mean(fg_bursts_occupancy_gatherer{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if less than 20 bursts are detected for either slow or fast gamma, then exclude that electrode
        if (length(length_gatherer_sg{i})<20) || (length(length_gatherer_fg{i})<20)
            mean_overlap_sg_fg_gatherer(i)=0;
            mean_overlap_sg_fg_gatherer_shuffle_trials(i)=0;
            mean_burst_occupancy_sg_gatherer(i)=0;
            mean_burst_occupancy_fg_gatherer(i)=0;
            power_gatherer_sg(i)=0;
            power_gatherer_fg(i)=0;
            continue;
        end
    end
    mean_overlap_sg_fg_gatherer=mean_overlap_sg_fg_gatherer(mean_overlap_sg_fg_gatherer~=0);
    mean_overlap_sg_fg_gatherer_shuffle_trials=mean_overlap_sg_fg_gatherer_shuffle_trials(mean_overlap_sg_fg_gatherer_shuffle_trials~=0);
    mean_burst_occupancy_sg_gatherer=mean_burst_occupancy_sg_gatherer(mean_burst_occupancy_sg_gatherer~=0);
    mean_burst_occupancy_fg_gatherer=mean_burst_occupancy_fg_gatherer(mean_burst_occupancy_fg_gatherer~=0);
    power_gatherer_sg=power_gatherer_sg(power_gatherer_sg~=0);
    power_gatherer_fg=power_gatherer_fg(power_gatherer_fg~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    power_gatherer_sg=10*log10(power_gatherer_sg);
    power_gatherer_fg=10*log10(power_gatherer_fg);
    subplot(plotHandles(Monkey_num,4));
    group=categorical({'Shuffled (trials)','Actual'});
    new_violin_swarm_plot_paired(mean_overlap_sg_fg_gatherer_shuffle_trials,mean_overlap_sg_fg_gatherer,0,1,group);
    ylabel('Mean Overlap (s)');
    if (Monkey_num==1)
        load('matched_power_idx_M1_fig_2.mat','matched_sg_indices','matched_fg_indices');
    elseif (Monkey_num==2)
        load('matched_power_idx_M2_fig_2.mat','matched_sg_indices','matched_fg_indices');
    end
    matched_sg_occupancy = mean_burst_occupancy_sg_gatherer(matched_sg_indices);
    matched_fg_occupancy = mean_burst_occupancy_fg_gatherer(matched_fg_indices);
    subplot(plotHandles(Monkey_num,3));
    mean_violin_swarm_plot(matched_sg_occupancy,matched_fg_occupancy);
    ylabel("Burst occupancy (s)");
end
  