% This script removes background noise from OD measurement and normalize
% across different measurements using PA14
% Last modified by Chen Liao, 04/24/2020

%% path to data
source_path = '../../data/growth_curve_original/';

%% files to read
filenames = {'010720GlycerolGrowthCurve_Hildi.csv';...
             '010720GlycerolGrowthCurve_Maurice.csv';...
             '011520GlycerolGrowthCurve_Hildi.csv';...
             '011520GlycerolGrowthCurve_Maurice.csv';...
             '012120GlycerolGrowthCurve_Hildi.csv';...
             '012120GlycerolGrowthCurve_Spark.csv';...
             '030420GlycerolGrowthCurve_Hildi.csv'};
nfiles = length(filenames);

%% strains for columns of 96-well plate
strains = {'Blank','F22031','F30658','F34365','F63912','F9670','M1608','H27930','H47921','H5708','PA14','Water';...
           'Blank','T38079','PA14','Water','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty';...
           'Blank','W70332','X9820','PA14','Water','Empty','Empty','Empty','Empty','Empty','Empty','Empty';...
           'Blank','M55212','M6057','M74707','T52373','T6313','T63266','W16407','W25637','W45909','PA14','Water';...
           'Blank','W91453','X78812','PA14','Water','Empty','Empty','Empty','Empty','Empty','Empty','Empty';...
           'Blank','F20590','F23197','F30658','F5677','M37351','S86968','T36994','W36662','W60856','PA14','Water';...
           'Blank','PA14','PAO1','PA7','Water','Empty','Empty','Empty','Empty','Empty','Empty','Empty'};
         
%% remove backrgound noise
index_of_background_column = [12, 4, 5, 12, 5, 12, 5];
growth_curve_bg_removed = cell(length(filenames), 3); % filename, time, OD
for i=1:length(filenames)
    
    % read raw growth curves
    fn_i = filenames{i};
    raw_data_i = readtable(strcat(source_path, fn_i));
    time_i = raw_data_i{:,2}/3600;
    plate_i = raw_data_i(:,[4:end]);
    
    % remove B4 in 010720GlycerolGrowthCurve_Hildi.csv
    if (i==1)
        plate_i(:,'B4') = [];
    end
    
    % remove growth curves at first and last rows (row A and H)
    plate_i(:,contains(plate_i.Properties.VariableNames,'A')) = [];
    plate_i(:,contains(plate_i.Properties.VariableNames,'H')) = [];
    
    % remove background noise
    strain_i = cell(size(plate_i,2),1);
    for j=1:size(plate_i,2)
        col_name = plate_i.Properties.VariableNames{j};
        letter_of_curr_col_name = col_name(1);
        background_well = char(strcat(letter_of_curr_col_name, string(index_of_background_column(i))));
        plate_i{:, col_name} = plate_i{:, col_name} - plate_i{:, background_well};
        strain_i{j} = strains{i,str2num(col_name(2:end))};
    end
    
    % create mean growth curve
    unique_strains_i = unique(strain_i);
    unique_strains_i(strcmp(unique_strains_i,'Water')) = [];
    unique_strains_i(strcmp(unique_strains_i,'Blank')) = [];
    unique_strains_i(strcmp(unique_strains_i,'Empty')) = [];
    
    variable_type    = cell(1, length(unique_strains_i));
    variable_type(:) = {'double'};
    plate_i_mean = table('Size',[size(plate_i,1),length(unique_strains_i)],'VariableType',variable_type,'VariableNames',unique_strains_i);
    for j=1:length(unique_strains_i)
        mean_curve_j = mean(plate_i{:,find(strcmp(strain_i, unique_strains_i(j)))}, 2);
        
        % if OD is non-positive at any point, replace it with the nearest positive value
        mean_curve_j(mean_curve_j<=0) = NaN;
        mean_curve_j = fillmissing(mean_curve_j, 'nearest');
        
        plate_i_mean{:,j} = mean_curve_j;
    end
    
    % store data
    growth_curve_bg_removed{i,1} = strrep(fn_i,'.csv','');
    growth_curve_bg_removed{i,2} = time_i;
    growth_curve_bg_removed{i,3} = plate_i_mean;
end

%% interpolation using new time points
timepoints = [0:0.1:48]'; % hour
ntps = length(timepoints);
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    
    variable_type    = cell(1, size(plate_i,2));
    variable_type(:) = {'double'};
    new_plate_i = table('Size',[ntps,size(plate_i,2)],'VariableType',variable_type,'VariableNames',plate_i.Properties.VariableNames);
    for j=1:size(plate_i,2)
        new_plate_i{:,j} = pchip(growth_curve_bg_removed{i,2}, plate_i{:,j}, timepoints);
    end
    
    growth_curve_bg_removed{i,2} = timepoints;
    growth_curve_bg_removed{i,3} = new_plate_i;
end

%% normalize growth curves based on PA14

% prepare table for linear mixed-effect model
stacked_plate_ids = cell(ntps*nfiles, 1);
stacked_timepoints = repmat(timepoints,nfiles,1);
stacked_od = zeros(ntps*nfiles, 1);
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    stacked_od((i-1)*ntps+1:i*ntps,1) = plate_i{:,'PA14'};
    stacked_plate_ids((i-1)*ntps+1:i*ntps,1) = {growth_curve_bg_removed{i,1}};
end
tblPA14 = table(stacked_plate_ids,stacked_timepoints,log10(stacked_od),'VariableNames',{'Plate','Time','log10OD'});
tblPA14.Time = categorical(tblPA14.Time);

% linear regressoin
opt_fit = fitlm(tblPA14, 'log10OD ~ Plate + Time', 'intercept', false, 'RobustOpts','logistic');

% plot adjusted PA14 growth curves vs. the originals
figure();
predicted_od = predict(opt_fit,tblPA14,'Prediction','Observation');
for i=1:nfiles
    subplot(2,4,i);
    hold on;
    plot(timepoints, 10.^(predicted_od((i-1)*ntps+1:i*ntps)), 'k--');
    plot(timepoints, stacked_od((i-1)*ntps+1:i*ntps), 'k-');
    ylabel('Predicted/Measured PA14');
    xlabel('Time (hour)');
    set(gca,'Xtick',[0:10:50]);
    set(gca,'Ytick',[0:0.5:1.5]);
    axis square;
    box on;
    legend({'adjusted';'original'}, 'Location', 'SouthEast');
    title(strrep(filenames{i},'.csv',''), 'interpreter', 'None');
end

% plot PA14 growth curves before adjustment
figure();
hold on;
my_colors = jet(nfiles);
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    plot(timepoints, plate_i{:,'PA14'}, '-', 'Color', my_colors(i,:));
end
ylabel('PA14 curves before ajustment');
xlabel('Time (hour)');
set(gca,'Xtick',[0:10:50]);
set(gca,'Ytick',[0:0.5:1.0]);
axis square;
box on;
legend({growth_curve_bg_removed{:,1}}, 'Interpreter', 'None', 'Location', 'NorthWest');

% adjust all growth curves using estimated coefficients
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    index = find(contains(opt_fit.CoefficientNames, strcat('Plate_',growth_curve_bg_removed{i,1})));
    if(~isempty(index))
        for j=1:size(plate_i,2)
            plate_i{:,j} = plate_i{:,j}/10^(opt_fit.Coefficients.Estimate(index));
        end
        growth_curve_bg_removed{i,3} = plate_i;
    end
end

% plot adjusted PA14 growth curves
figure();
hold on;
my_colors = jet(nfiles);
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    plot(timepoints, plate_i{:,'PA14'}, '-', 'Color', my_colors(i,:));
end
ylabel('PA14 curves after ajustment');
xlabel('Time (hour)');
set(gca,'Xtick',[0:10:50]);
set(gca,'Ytick',[0:0.5:1.0]);
axis square;
box on;
legend({growth_curve_bg_removed{:,1}}, 'Interpreter', 'None', 'Location', 'NorthWest');

%% Compile a single table of growth curves including all strains
matgcall_w_duplicates = [];
strain_names_of_matgcall = [];
for i=1:nfiles
    plate_i = growth_curve_bg_removed{i,3};
    matgcall_w_duplicates = [matgcall_w_duplicates, plate_i{:,:}];
    strain_names_of_matgcall = [strain_names_of_matgcall, plate_i.Properties.VariableNames];
end

unique_strain_names = unique(strain_names_of_matgcall);
variable_type    = cell(1, length(unique_strain_names));
variable_type(:) = {'double'};
tbl_gc_all = table('Size',[ntps,length(unique_strain_names)],'VariableType',variable_type,'VariableNames',unique_strain_names);

for i=1:length(unique_strain_names)
    tbl_gc_all{:,unique_strain_names{i}} = mean(matgcall_w_duplicates(:,find(strcmp(strain_names_of_matgcall, unique_strain_names(i)))), 2);
end
tbl_gc_all.Properties.RowNames = cellstr(num2str(timepoints));

%% Write to file
writetable(tbl_gc_all, 'normalized_mean_growth_curve_PA_glycerol.csv', 'Delimiter', ',', 'WriteRowNames',true);