function [sg_bursts_union_length, fg_bursts_union_length, sg_fg_overlap_length, ...
frac_sg_overlapping, frac_fg_overlapping] = get_overlap_bursts(start_point_temp_sg,...
end_point_temp_sg, start_point_temp_fg, end_point_temp_fg)
    %%%%%%%%%%%%%%%%%%%%%%%%%% Slow gamma %%%%%%%%%%%%%%%%%%%%%%
    sg_bursts_interval = [start_point_temp_sg(:) end_point_temp_sg(:)];
    if ~isempty(sg_bursts_interval)
        sg_bursts_interval = sortrows(sg_bursts_interval,1);
        sg_bursts_interval_merge = sg_bursts_interval(1,:);
        for idx = 2:size(sg_bursts_interval,1)
            if sg_bursts_interval(idx,1) <= sg_bursts_interval_merge(end,2)
                sg_bursts_interval_merge(end,2) = max(sg_bursts_interval_merge(end,2), ...
                                       sg_bursts_interval(idx,2));
            else
                sg_bursts_interval_merge = [sg_bursts_interval_merge; sg_bursts_interval(idx,:)];
            end
        end
        sg_bursts_union_length = sum(sg_bursts_interval_merge(:,2)-sg_bursts_interval_merge(:,1));
    else
        sg_bursts_interval_merge=[];
        sg_bursts_union_length = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%% Fast gamma %%%%%%%%%%%%%%%%%%%%%%
    fg_bursts_interval = [start_point_temp_fg(:) end_point_temp_fg(:)];
    if ~isempty(fg_bursts_interval)
        fg_bursts_interval = sortrows(fg_bursts_interval,1);
        fg_bursts_interval_merge = fg_bursts_interval(1,:);
        for idx = 2:size(fg_bursts_interval,1)
            if fg_bursts_interval(idx,1) <= fg_bursts_interval_merge(end,2)
               fg_bursts_interval_merge(end,2) = max(fg_bursts_interval_merge(end,2), ...
                                       fg_bursts_interval(idx,2));
            else
                fg_bursts_interval_merge = [fg_bursts_interval_merge;...
                    fg_bursts_interval(idx,:)];
            end
        end
        fg_bursts_union_length = sum(fg_bursts_interval_merge(:,2)-fg_bursts_interval_merge(:,1));
    else
        fg_bursts_interval_merge=[];
        fg_bursts_union_length = 0;
    end
    %%%%%%%% Slow gamma and Fast gamma overlap duration %%%%%%%%%%%%%
    if isempty(sg_bursts_interval_merge) || isempty(fg_bursts_interval_merge)
        sg_fg_overlap_length = 0;
    else
        overlap_intervals = [];
        for idx_sg = 1:size(sg_bursts_interval_merge,1)
            for idx_fg = 1:size(fg_bursts_interval_merge,1)
                start_temp = max(sg_bursts_interval_merge(idx_sg,1), fg_bursts_interval_merge(idx_fg,1));
                end_temp = min(sg_bursts_interval_merge(idx_sg,2), fg_bursts_interval_merge(idx_fg,2));
                if end_temp > start_temp
                    overlap_intervals = [overlap_intervals; start_temp end_temp];
                end
            end
        end
        if ~isempty(overlap_intervals)
            overlap_intervals = sortrows(overlap_intervals,1);
            overlap_merged = overlap_intervals(1,:);
            for idx = 2:size(overlap_intervals,1)
                if overlap_intervals(idx,1) <= overlap_merged(end,2)
                    overlap_merged(end,2) = max(overlap_merged(end,2), ...
                                               overlap_intervals(idx,2));
                else
                    overlap_merged = [overlap_merged; overlap_intervals(idx,:)];
                end
            end
            sg_fg_overlap_length = sum(overlap_merged(:,2)-overlap_merged(:,1));
        else
            sg_fg_overlap_length = 0;
        end
    end
    %%%%%%% Fractions %%%%%%%%%%%%
    if sg_bursts_union_length == 0
        frac_sg_overlapping = [];
    else
        frac_sg_overlapping = sg_fg_overlap_length / sg_bursts_union_length;
    end
    if fg_bursts_union_length == 0
        frac_fg_overlapping = [];
    else
        frac_fg_overlapping = sg_fg_overlap_length / fg_bursts_union_length;
    end