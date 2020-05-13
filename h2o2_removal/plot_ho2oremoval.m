% This script plots OD and H2O2 removal of Pseudomonas strains
% Last updated by Chen Liao, 05/08/2020

%% Read data
rawData_OD = xlsread('../data/h2o2/181130_glycerol_Amplex.xlsx','OD600');
rawData_H2O2 = xlsread('../data/h2o2/181130_glycerol_Amplex.xlsx','AmplexEM');
time_48 = rawData_OD(:,2)/3600;  % in unit of hour
rawData_OD = rawData_OD(:,[4:end]);
rawData_H2O2 = rawData_H2O2(:,[4:end]);

%% strains
strain_id = {'Control','PA14','W25637','W36662','M1608','W91453','F22031','F9670','M37351','T38079','X9820','Control'};
mycolors = {'';'r';'r';'b';'b';'r';'r';'r';'g';'r';'r';''};
figure();

%% plot OD600
subplot(2,1,1);
hold on;
for i=2:11 % strain 
    % 2-4, w/ dye; 5-7, w/o dye
    plot(time_48,mean(rawData_OD(:,[(i-1)*8+2:(i-1)*8+4]),2),'k-','Color',mycolors{i},'LineWidth',1.5);
end
%axis square;
box on;
axis([0,48,0,1.5]);
set(gca,'XTick',[0:12:48]);
set(gca,'YTick',[0:0.5:1.5]);
xlabel('Time (h)');
ylabel('OD');

%% plot H2O2 removed using Amplex EM
subplot(2,1,2);
background = mean(rawData_H2O2(:,[2:4]),2);
for i=2:11 % strain
    hold on;
    % 2-4, w/ dye; 5-7, w/o dye
    plot(time_48,background-mean(rawData_H2O2(:,[(i-1)*8+2:(i-1)*8+4]),2),'k-','Color',mycolors{i},'LineWidth',1.5);
end
%axis square;
box on;
axis([0,48,-1000,1800]);
set(gca,'XTick',[0:12:48]);
set(gca,'YTick',[0:600:1800]);
xlabel('Time (h)');
ylabel('Absorbed flurescence');