% This script divides each growth curve into phase I, II and III
% last modified by Chen Liao, April 25, 2020

%% read data
tblgc = readtable('../data/growth_curve_normalized/normalized_mean_growth_curve_PA_succinate.csv');
tblgc.Properties.VariableNames{1} = 'Time';
timepoints = tblgc{:,1};
ntps = length(timepoints);
tblgc.Time = [];

%% identify growth phases
phase_start_index = zeros(size(tblgc,2),3);
phase_start_time = zeros(size(tblgc,2),3);
phase_end_time = zeros(size(tblgc,2),3);
for i=1:size(tblgc,2)
    
    % calculate first-order derivative
    yobs = tblgc{:,i}; % observed data
    yobs_prime = ppval(fnder(spline(timepoints, yobs),1), timepoints); % first order derivative
    yobs_prime_smoothed = sgolayfilt(yobs_prime, 3, 51); % smoothing first-order derivative
    
    % find growth phases
    index_pos_val = find(yobs_prime_smoothed>0.001);
    if (~isempty(index_pos_val))
        
        % Phase I & II: find the longest interval that first derivative is above 0
        max_interval_length = 0;
        curr_left_end = index_pos_val(1);
        curr_right_end = index_pos_val(1);
        phase1_start_index = curr_left_end;
        phase3_start_index = curr_right_end;
        for j=2:length(index_pos_val)
            % some temporary drops are ignored
            max_gap_length = 10;
            if ((index_pos_val(j)-index_pos_val(j-1))<max_gap_length)
                curr_right_end = index_pos_val(j);
                if (j==length(index_pos_val) & curr_right_end-curr_left_end > max_interval_length)
                    phase1_start_index = curr_left_end;
                    phase3_start_index = curr_right_end;
                end
            else
                if (curr_right_end-curr_left_end > max_interval_length)
                    
                    % special treatment
                   
                    if (i==10)
                        % H5708
                        if (timepoints(curr_right_end)>20)
                            break
                        end
                    end
                    if (i==17)
                        % PA7
                        if (timepoints(curr_right_end)>20)
                            break
                        end
                    end
                    if (i==23)
                        % T63266
                        if (timepoints(curr_right_end)>20)
                            break
                        end
                    end
                    
                    phase1_start_index = curr_left_end;
                    phase3_start_index = curr_right_end;
                    max_interval_length = phase3_start_index - phase1_start_index;
                end
                curr_left_end = index_pos_val(j);
                curr_right_end = index_pos_val(j);
            end
        end

        phase_start_index(i,1) = phase1_start_index;
        phase_start_index(i,3) = phase3_start_index;
        phase_start_time(i,1) = timepoints(phase1_start_index);
        phase_start_time(i,3) = timepoints(phase3_start_index);
        phase_end_time(i,2) = timepoints(phase3_start_index);
        
        % Phase II: if multiple peaks exist within the interval, choose the
        % highest peak as the entry point for Phase II
        [pks, locs] = findpeaks(yobs_prime_smoothed);
        pks = pks(find(locs >= phase1_start_index & locs <= phase3_start_index));
        locs = locs(find(locs >= phase1_start_index & locs <= phase3_start_index));
        if (~isempty(locs))
            [~,max_peak_index] = max(pks);
            phase_start_index(i,2) = locs(max_peak_index);
            phase_start_time(i,2) = timepoints(locs(max_peak_index));
        else
            phase_start_index(i,2) = NaN;
            phase_start_time(i,2) = NaN;
        end
        
        % find end time of phase III
        
    else
        fprintf('%s: no growth at all.\n', tblgc.Properties.VariableNames{i});
        phase_start_index(i,:) = [NaN, NaN, NaN];
        phase_start_time(i,:) = [NaN, NaN, NaN];
    end
end

%% plot growth phases
figure();
for i=1:size(tblgc,2)
    subplot(5,7,i);
    hold on;
    
    yobs = tblgc{:,i};  
    plot(timepoints, yobs, 'k-', 'LineWidth', 1);
    
    if (sum(isnan(phase_start_index(i,:)))==0)
        phase1_start_index = phase_start_index(i,1);
        phase2_start_index = phase_start_index(i,2);
        phase3_start_index = phase_start_index(i,3);
        patchline(timepoints(phase1_start_index:phase2_start_index), yobs(phase1_start_index:phase2_start_index),...
            'linestyle', '-', 'edgecolor', 'r', 'linewidth', 4, 'edgealpha', 0.2);
        patchline(timepoints(phase2_start_index:phase3_start_index), yobs(phase2_start_index:phase3_start_index),...
            'linestyle', '-', 'edgecolor', 'b', 'linewidth', 4, 'edgealpha', 0.2);
        patchline(timepoints(phase3_start_index:end), yobs(phase3_start_index:end),...
            'linestyle', '-', 'edgecolor', 'g', 'linewidth', 4, 'edgealpha', 0.2);
    end
    
    xlim([0,48]);
    set(gca,'XTick',[0,12,24,36,48]);
    ylim([0,0.4]);
    set(gca,'YTick',[0,0.2,0.4]);
    axis square;
    box on;
    xlabel('Time (hour)');
    ylabel('OD');
    title(tblgc.Properties.VariableNames{i});
end

%% save to file
tbl_phase_start_time = array2table(phase_start_time);
tbl_phase_start_time.Properties.RowNames = tblgc.Properties.VariableNames;
tbl_phase_start_time.Properties.VariableNames = {'Phase1';'Phase2';'Phase3'};
writetable(tbl_phase_start_time, 'PA_succinate_growth_phase_start_time.csv', 'Delimiter', ',', 'WriteRowNames', true);