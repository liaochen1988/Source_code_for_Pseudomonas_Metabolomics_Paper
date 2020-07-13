% This script plots OD and H2O2 removal of Pseudomonas strains
% Last updated by Chen Liao, 05/08/2020

%% Read data
rawData_OD = xlsread('../data/h2o2/181130_glycerol_Amplex.xlsx','OD600');
rawData_EM = xlsread('../data/h2o2/181130_glycerol_Amplex.xlsx','AmplexEM');
time_48 = rawData_OD(:,2)/3600;  % in unit of hour
rawData_OD = rawData_OD(:,[4:end]);
rawData_EM = rawData_EM(:,[4:end]);

%% strains
strain_id = {'Control','PA14','W25637','W36662','M1608','W91453','F22031','F9670','M37351','T38079','X9820','Control'};
mycolors = {'';'r';'r';'b';'b';'r';'r';'r';'g';'r';'r';''};
lw=1.5;
figure();

%% calculate OD, accumulated removal of H2O2, and per OD removal rate
order = 3;
framelen = 21;

OD_corrected = zeros(length(timeOD),12,3);
derivOD_corrected = zeros(length(timeOD),12,3);
gRate_corrected = zeros(length(timeOD),12,3);
EM_corrected = zeros(length(timeEM),12,3);
derivEM_corrected = zeros(length(timeEM),12,3);
rmRate_corrected = zeros(length(timeEM),12,3);

for i=2:11 % strain
    
    % OD600
    for k=1:3
        OD_corrected(:,i,k) = rawDataOD(:,(i-1)*8+k+1)-rawDataOD(:,k+1);
    end
    
    % OD derivative
    for k=1:3
        pp=spline(timeOD,sgolayfilt(squeeze(OD_corrected(:,i,k)),order,framelen));
        p_der=fnder(pp,1);
        derivOD_corrected(:,i,k) = sgolayfilt(ppval(p_der,timeOD),order,framelen);
    end
    
    % Growth rate
    for k=1:3
        gRate_corrected(:,i,k) = sgolayfilt(derivOD_corrected(:,i,k)./OD_corrected(:,i,k),order,framelen);
    end
    
    % EM
    for k=1:3
        EM_corrected(:,i,k) = rawDataEM(:,k+1)-rawDataEM(:,(i-1)*8+k+1);
    end
    
    % EM derivative (d EM/dt )
    for k=1:3
        pp=spline(timeEM,sgolayfilt(squeeze(EM_corrected(:,i,k)),order,framelen));
        p_der=fnder(pp,1);
        derivEM_corrected(:,i,k) = sgolayfilt(ppval(p_der,timeOD),order,framelen);
    end
    
    % H2O2 removal rate (d EM/dt/OD)
    for k=1:3
        rmRate_corrected(:,i,k) = sgolayfilt(derivEM_corrected(:,i,k)./OD_corrected(:,i,k),order,framelen);
    end
    
end

%% plot OD
subplot(3,1,1);
hold on;
for i=2:11
    if (i<11)
        OD_mean = mean(squeeze(OD_corrected(:,i,:)),2);
        OD_std = std(squeeze(OD_corrected(:,i,:)),0,2);
    else
        OD_mean = mean(squeeze(OD_corrected(:,i,[2,3])),2);
        OD_std = std(squeeze(OD_corrected(:,i,[2,3])),0,2);
    end
    OD_ub = OD_mean+OD_std;
    OD_lb = OD_mean-OD_std;
    plot(timeOD,OD_mean,'k-','Color',mycolors{i},'LineWidth',lw);
    h=fill([timeOD;flip(timeOD)],[OD_lb;flip(OD_ub)],mycolors{i});
    set(h,'facealpha',.1);
end
box on;
axis([0,48,0,1.2]);
set(gca,'XTick',[0:12:48]);
set(gca,'YTick',[0:0.4:1.2]);
xlabel('Time (h)');
ylabel('OD');

%% plot accumulated H2O2 removal
subplot(3,1,2);
hold on;
for i=2:11
    if (i<11)
        EM_mean = mean(squeeze(EM_corrected(:,i,:)),2);
        EM_std = std(squeeze(EM_corrected(:,i,:)),0,2);
    else
        EM_mean = mean(squeeze(EM_corrected(:,i,[2,3])),2);
        EM_std = std(squeeze(EM_corrected(:,i,[2,3])),0,2);
    end
    EM_ub = EM_mean+EM_std;
    EM_lb = EM_mean-EM_std;
    plot(timeOD,EM_mean,'k-','Color',mycolors{i},'LineWidth',lw);
    h=fill([timeOD;flip(timeOD)],[EM_lb;flip(EM_ub)],mycolors{i});
    set(h,'facealpha',.1);
end
box on;
axis([0,48,-600,1800]);
set(gca,'XTick',[0:12:48]);
set(gca,'YTick',[-600:600:1800]);
xlabel('Time (h)');
ylabel('H_2O_2 removal');

%% plot single cell detoxification rate
subplot(3,1,3);
hold on;
for i=2:11
    if (i<11)
        rmRate_mean = mean(squeeze(rmRate_corrected(:,i,:)),2);
        rmRate_std = std(squeeze(rmRate_corrected(:,i,:)),0,2);
    else
        rmRate_mean = mean(squeeze(rmRate_corrected(:,i,[2,3])),2);
        rmRate_std = std(squeeze(rmRate_corrected(:,i,[2,3])),0,2);
    end
    rmRate_ub = rmRate_mean+rmRate_std;
    rmRate_lb = rmRate_mean-rmRate_std;
    plot(timeOD,rmRate_mean,'k-','Color',mycolors{i},'LineWidth',lw);
    h=fill([timeOD;flip(timeOD)],[rmRate_lb;flip(rmRate_ub)],mycolors{i});
    set(h,'facealpha',.1);
end
box on;
axis([0,48,-3200,3200]);
set(gca,'XTick',[0:12:48]);
set(gca,'YTick',[-2400:1200:2400]);
xlabel('Time (h)');
ylabel('H_2O_2 removal rate per OD');
