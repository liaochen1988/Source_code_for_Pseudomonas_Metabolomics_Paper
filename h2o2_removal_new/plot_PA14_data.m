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
 
%% remove backrgound h2o2
figure();
hold on;
my_colors = {'r';'b';'g';'k'};

index_of_background_column = [1, 1, 1, 1];
for i=1:length(filenames)
    
    % read raw growth curves
    fn_i = filenames{i};
    raw_data_em_i = readtable(strcat(source_path, fn_i),'Sheet','AmplexEM');
    time_em_i = raw_data_em_i{:,2}/3600;
    em_i = raw_data_em_i(:,[4:end]);
        
    % only keep row B,C,D
    rows_to_remove = {'E';'F';'G';'H'};
    for j=1:length(rows_to_remove)
        em_i(:,contains(em_i.Properties.VariableNames, rows_to_remove{j})) = [];
    end
    
    % remove background noise
    for j=size(em_i,2):-1:1
        col_name = em_i.Properties.VariableNames{j};
        letter_of_curr_col_name = col_name(1);
        background_well = char(strcat(letter_of_curr_col_name, string(index_of_background_column(i))));
        em_i{:, col_name} = em_i{:, background_well} - em_i{:, col_name};
    end
      
    plot(time_em_i, em_i{:,'A2'}, '-', 'Color', my_colors{i});
    plot(time_em_i, em_i{:,'B2'}, '-', 'Color', my_colors{i});
    plot(time_em_i, em_i{:,'C2'}, '-', 'Color', my_colors{i});
    plot(time_em_i, em_i{:,'D2'}, '-', 'Color', my_colors{i});
end

