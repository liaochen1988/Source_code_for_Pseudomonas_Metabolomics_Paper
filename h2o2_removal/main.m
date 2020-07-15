% This script plots OD and H2O2 removal of Pseudomonas strains, and run a
% mixed effect linear model to test the effect of rhamnolipid production on
% H2O2 removal capability
%
% Last updated by Chen Liao, 07/14/2020

%% files to read
source_path = '../data/h2o2/';
filenames = {'062420H2O2Removal_Hildi.xls';...
             '063020H2O2Production_Hildi.xls';...
             '063020H2O2Production_Maurice.xls';...
             '070620H2O2Production_Hildi.xls'};
nfiles = length(filenames);
 
%% strains for columns of 96-well plate
strains = {'PBS','PA14','F30658','F34365','F63912','H27930','H47921','H5708','M6075','M74707','S86968','H2O2';...
           'PBS','PA14','T63266','W16407','W45909','W60856','W70332','X78812','H2O2','Empty','Empty','Empty';...
           'PBS','PA14','F23197','F5677','M55212','PA7','PAO1','T52373','T6313','H2O2','Empty','Empty';...
           'PBS','PA14','F22031','F9670','M1608','M37351','T38079','W25637','W36662','W91453','X9820','H2O2'};
         
%% find pa14 average
timepoints = [0:1/6:24]'; % hour, 10 min interval
PA14_ave = zeros(length(timepoints), length(filenames));
for i=1:length(filenames)
    
    % read raw growth curves
    fn_i = filenames{i};
    raw_data_em_i = readtable(strcat(source_path, fn_i),'Sheet','AmplexEM');
    time_em_i = raw_data_em_i{:,2}/3600;
    em_i = raw_data_em_i(:,[4:end]);
    PA14_ave(:,i) = pchip(time_em_i, mean(em_i{:, {'B1';'C1';'D1'}}-em_i{:,{'B2';'C2';'D2'}},2), timepoints);
end

%% normalize h2o2 removal relative to pa14
index_of_background_column = [1, 1, 1, 1];
normalized_data = cell(length(filenames), 3); % filename, OD, EM
for i=1:length(filenames)
    
    % read raw growth curves
    fn_i = filenames{i};
    raw_data_od_i = readtable(strcat(source_path, fn_i),'Sheet','OD');
    raw_data_em_i = readtable(strcat(source_path, fn_i),'Sheet','AmplexEM');

    time_od_i = raw_data_od_i{:,2}/3600;
    od_i = raw_data_od_i(:,[4:end]);
    time_em_i = raw_data_em_i{:,2}/3600;
    em_i = raw_data_em_i(:,[4:end]);
        
    % only keep row B,C,D
    rows_to_remove = {'A';'E';'F';'G';'H'};
    for j=1:length(rows_to_remove)
        od_i(:,contains(od_i.Properties.VariableNames, rows_to_remove{j})) = [];
        em_i(:,contains(em_i.Properties.VariableNames, rows_to_remove{j})) = [];
    end
    
    % remove background noise
    od_i_new = array2table(zeros(length(timepoints), size(od_i,2)), 'VariableNames', od_i.Properties.VariableNames);
    em_i_new = array2table(zeros(length(timepoints), size(em_i,2)), 'VariableNames', od_i.Properties.VariableNames);
    strain_i = cell(size(od_i,2),1);
    for j=size(od_i,2):-1:1
        col_name = od_i.Properties.VariableNames{j};
        letter_of_curr_col_name = col_name(1);
        background_well = char(strcat(letter_of_curr_col_name, string(index_of_background_column(i))));
        od_corrected = od_i{:, col_name} - od_i{:, background_well};
        od_corrected(od_corrected<0) = NaN; % if OD is non-positive at any point, replace it with the nearest positive value
        od_i_new{:, col_name} = pchip(time_od_i, fillmissing(od_corrected, 'nearest'), timepoints);
        em_i_new{:, col_name} = pchip(time_em_i, em_i{:, background_well} - em_i{:, col_name}, timepoints)./PA14_ave(:,i);
        
        switch col_name(1)
            case 'B'
                replicate='_R1';
            case 'C'
                replicate='_R2';
            case 'D'
                replicate='_R3';
            otherwise
                error('unknown column name');
        end
        fn_i_split = split(fn_i,'H2O2');
        date = fn_i_split(1);
        fn_i_split = split(fn_i_split(2), '_');
        machine = strrep(fn_i_split(2),'.xls','');
        strain_i{j} = strcat(strains{i,str2num(col_name(2:end))}, replicate, '_', date{1}, '_', machine{1});
    end
    od_i_new.Properties.VariableNames = strain_i;
    od_i_new(:,contains(od_i_new.Properties.VariableNames,'PBS')) = [];
    od_i_new(:,contains(od_i_new.Properties.VariableNames,'Empty')) = [];
    od_i_new(:,contains(od_i_new.Properties.VariableNames,'H2O2')) = [];

    em_i_new.Properties.VariableNames = strain_i;
    em_i_new(:,contains(em_i_new.Properties.VariableNames,'PBS')) = [];
    em_i_new(:,contains(em_i_new.Properties.VariableNames,'Empty')) = [];
    em_i_new(:,contains(em_i_new.Properties.VariableNames,'H2O2')) = [];
    
    % store data
    normalized_data{i,1} = strrep(fn_i,'.csv','');
    normalized_data{i,2} = od_i_new;
    normalized_data{i,3} = em_i_new;
end

%% Compile a single table of all strain/replicates
normalized_od_all = [];
normalized_em_all = [];
for i=1:nfiles
    od_i = normalized_data{i,2};
    em_i = normalized_data{i,3};
    
    if (i==1)
        normalized_od_all = [normalized_od_all, od_i];
        normalized_em_all = [normalized_em_all, em_i];
    else
        normalized_od_all = [normalized_od_all, od_i(:,~contains(od_i.Properties.VariableNames,'PA14'))];
        normalized_em_all = [normalized_em_all, em_i(:,~contains(em_i.Properties.VariableNames,'PA14'))];
    end
end

%% read rhamnolipid production phenotype
RL_phenotype = readtable('../data/rhamnolipids/rhamnMat.xlsx');

%% plot experimental data
order = 3;
framelen = 11;
mycolors = {'b';'g';'r'}; % blue: non-producer, green: mild producer, red: strong producer

figure();
rmRate = normalized_em_all;
for j=1:size(normalized_od_all,2) 
    varname = normalized_od_all.Properties.VariableNames{j};
    str_split = split(varname, '_');
    strain_name = str_split{1};
    if (strcmp(strain_name, 'M6075'))
        continue
    end
    rhl_level = RL_phenotype{strcmp(RL_phenotype.strain, strain_name), 'rhamn3cats'}+1;
    
    % OD
    subplot(3,1,1);
    hold on;
    plot(timepoints,normalized_od_all{:,j},'k-','Color',mycolors{rhl_level},'LineWidth',1);
    box on;
    axis([0,24,0,1]);
    set(gca,'XTick',[0:6:24]);
    set(gca,'YTick',[0:0.2:1]);
    xlabel('Time (h)');
    ylabel('OD');
    
    % H2O2 removal
    subplot(3,1,2);
    hold on;
    plot(timepoints,normalized_em_all{:,j},'k-','Color',mycolors{rhl_level},'LineWidth',1);
    box on;
    axis([0,24,-0.5,1.5]);
    set(gca,'XTick',[0:6:24]);
    set(gca,'YTick',[-0.5:0.5:1.5]);
    xlabel('Time (h)');
    ylabel('H_2O_2 removal');

    % H2O2 removal rate
    pp=spline(timepoints,sgolayfilt(normalized_em_all{:,j},order,framelen));
    p_der=fnder(pp,1);
    derivEM = sgolayfilt(ppval(p_der,timepoints),order,framelen);
    rmRate{:,j} = sgolayfilt(derivEM./normalized_od_all{:,j},order,framelen);
    
    subplot(3,1,3);
    hold on;
    plot(timepoints,rmRate{:,j},'k-','Color',mycolors{rhl_level},'LineWidth',1);
    box on;
    axis([0,24,-40,20]);
    set(gca,'XTick',[0:6:24]);
    set(gca,'YTick',[-40:20:20]);
    xlabel('Time (h)');
    ylabel('H_2O_2 removal rate per OD');
end

%% identify effect of RL production on h2o2 removal
tblreg = array2table(zeros(size(rmRate,1)*size(rmRate,2),4),'VariableNames',{'RemovalRate','Time','RL','CurveID'});
for j=1:size(rmRate,2)
    varname = normalized_od_all.Properties.VariableNames{j};
    str_split = split(varname, '_');
    strain_name = str_split{1};
    rhl_level = RL_phenotype{strcmp(RL_phenotype.strain, strain_name), 'rhamn3cats'};

    tblreg{((j-1)*size(rmRate,1)+1):j*size(rmRate,1),1} = rmRate{:,j};
    tblreg{((j-1)*size(rmRate,1)+1):j*size(rmRate,1),2} = timepoints;
    tblreg{((j-1)*size(rmRate,1)+1):j*size(rmRate,1),3} = rhl_level*ones(size(rmRate,1),1);
    tblreg{((j-1)*size(rmRate,1)+1):j*size(rmRate,1),4} = j*ones(size(rmRate,1),1);
end
tblreg.Time = categorical(tblreg.Time);
tblreg.RL = categorical(tblreg.RL);
tblreg.CurveID = categorical(tblreg.CurveID);

% RL_1, estimate 1.2479, p-value 3.08e-4
% RL_2, estimate 1.4109, p-value 8.49e-7
lme = fitlme(tblreg,'RemovalRate~Time+RL+(1|CurveID)');