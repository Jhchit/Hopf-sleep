%%%% code for Fig. 5

clear all;
close all;

%% data loading
filename = 'Hopf-Sleep\analysis\data for analysis';
load([filename,'\ts.mat']); % file of bold signal 
G_A = xlsread([filename,'\gender_age.xlsx']); % file of gender and age
AGE = G_A(:,3);
GENDER = G_A(:,2);

%% time-lagged FC calculation
flag = 2;% 1 --- calculate FC without taking absolute value;
         % 2 --- calculate FC then take absolute value;
exclude = 1; % 0 --- not exclude outlier
             % 1 --- exclude outlier
excluding_method.method = 'std' ; % is --- isoutlier, std --- standard deviation
excluding_method.para = 1.5;

%  SLO
 % morning-SLO
 data = ts.A.m;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3); % hypothalamus
     % ACC-basal forebrain
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_m_A(i,1) = I-21;
     
     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_m_A(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_m_A(i,1) = I-21;
 end

 if exclude == 1
     peak_timing_A_BF_m_A = exclude_outlier(peak_timing_A_BF_m_A,excluding_method);
     peak_timing_A_TH_m_A = exclude_outlier(peak_timing_A_TH_m_A,excluding_method);
     peak_timing_BF_TH_m_A = exclude_outlier(peak_timing_BF_TH_m_A,excluding_method);
 end




 % night-SLO
 data = ts.A.n;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3);%hypothalamus
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_n_A(i,1) = I-21;

     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_n_A(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_n_A(i,1) = I-21;
 end

 if exclude == 1
     peak_timing_A_BF_n_A = exclude_outlier(peak_timing_A_BF_n_A,excluding_method);
     peak_timing_A_TH_n_A = exclude_outlier(peak_timing_A_TH_n_A,excluding_method);
     peak_timing_BF_TH_n_A = exclude_outlier(peak_timing_BF_TH_n_A,excluding_method);
 end

%  SLW
 % morning-SLW
 data = ts.B.m;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3); % hypothalamus
     % ACC-Basal forebrain
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_m_B(i,1) = I-21;

     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_m_B(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_m_B(i,1) = I-21;
 end

 if exclude == 1
     peak_timing_A_BF_m_B = exclude_outlier(peak_timing_A_BF_m_B,excluding_method);
     peak_timing_A_TH_m_B = exclude_outlier(peak_timing_A_TH_m_B,excluding_method);
     peak_timing_BF_TH_m_B = exclude_outlier(peak_timing_BF_TH_m_B,excluding_method);
 end


 % night-SLW
 data = ts.B.n;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3); % hypothalamus
     % ACC-Basal forebrain
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_n_B(i,1) = I-21;

     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_n_B(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_n_B(i,1) = I-21;
 end

 if exclude == 1
     peak_timing_A_BF_n_B = exclude_outlier(peak_timing_A_BF_n_B,excluding_method);
     peak_timing_A_TH_n_B = exclude_outlier(peak_timing_A_TH_n_B,excluding_method);
     peak_timing_BF_TH_n_B = exclude_outlier(peak_timing_BF_TH_n_B,excluding_method);
 end

 %  HC
 % morning-HC
 data = ts.C.m;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3); % hypothalamus
     % ACC-Basal forebrain
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_m_C(i,1) = I-21;

     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_m_C(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_m_C(i,1) = I-21;
 end

 if exclude == 1
     peak_timing_A_BF_m_C = exclude_outlier(peak_timing_A_BF_m_C,excluding_method);
     peak_timing_A_TH_m_C = exclude_outlier(peak_timing_A_TH_m_C,excluding_method);
     peak_timing_BF_TH_m_C = exclude_outlier(peak_timing_BF_TH_m_C,excluding_method);
 end


 % night-HC
 data = ts.C.n;
 for i = 1:length(data)
     sig = data{i};
     ACC = sig(:,1); % acc
     BF = sig(:,2); % basal forebrain
     TH = sig(:,3); % hypothalamus
     % ACC-Basal forebrain
     [r_lag] = Pearson_time_lage(ACC,BF,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_BF_n_C(i,1) = I-21;

     % ACC-Hypothalamus
     [r_lag] = Pearson_time_lage(ACC,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_A_TH_n_C(i,1) = I-21;

     % Basal forebrain-Hypothalamus
     [r_lag] = Pearson_time_lage(BF,TH,20);
     if flag == 2
         r_lag = abs(r_lag);
     end
     [M,I] = max(r_lag);
     peak_timing_BF_TH_n_C(i,1) = I-21;
 end
 if exclude == 1
     peak_timing_A_BF_n_C = exclude_outlier(peak_timing_A_BF_n_C,excluding_method);
     peak_timing_A_TH_n_C = exclude_outlier(peak_timing_A_TH_n_C,excluding_method);
     peak_timing_BF_TH_n_C = exclude_outlier(peak_timing_BF_TH_n_C,excluding_method);
 end

 figure(7);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,50,600,690]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .65]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
RGB1 = [0.80392,0.78824,0.78824]; % light grey[0.28235,0.23922,0.5451]; % highly infiltrated regions[0.86667,0.62745,0.86667]; % light purple[0.80392,0.78824,0.78824]; % light grey[0.5451,0.53725,0.53725]; % dark grey[0.80392,0.78824,0.78824]; % light grey for low-grade%[0.00, 0.45, 0.74];%[0.00,0.45,0.74];%蓝色  [217,205,204]/264;
RGB2 = [0.5451,0.53725,0.53725]; % dark grey[0.8549,0.64706,0.12549]; % lowly infiltrated regions[0.84706,0.74902,0.84706]; % light red[1,0.89412,0.7098];% healthy%[0.85, 0.33, 0.10];%[0.85,0.33,0.10];%红色  [122,150,171]/264;
RGB3 = [0.91373,0.58824,0.47843]; % surrounding regions[1,0.71373,0.75686]; % light pink [0.91373,0.58824,0.47843]; % surrounding regions[0.5451,0.53725,0.53725]; % dark grey high-grade [0.82353,0.41176,0.11765]; % surrounding regions[191,141,122]/264;%[0.9290 0.6940 0.1250];
RGB4 = [1,0.89412,0.7098]; % light red
RGB = [RGB1;RGB2;RGB3;RGB4];

%% plot
figure(1);
subplot(321);
% A = struct('BF_ACC',-peak_timing_A_BF_m_C,'BF_TH',peak_timing_BF_TH_m_C,...
%     'TH_ACC',-peak_timing_A_TH_m_C,'TH_BF',-peak_timing_BF_TH_m_C);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_m_C;peak_timing_BF_TH_m_C;...
    -peak_timing_A_TH_m_C;-peak_timing_BF_TH_m_C};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('Healthy, Morning');

subplot(322);
% A = struct('BF_ACC',-peak_timing_A_BF_n_C,'BF_TH',peak_timing_BF_TH_n_C,...
%     'TH_ACC',-peak_timing_A_TH_n_C,'TH_BF',-peak_timing_BF_TH_n_C);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_n_C;peak_timing_BF_TH_n_C;...
    -peak_timing_A_TH_n_C;-peak_timing_BF_TH_n_C};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('Healthy, Evening');

subplot(323);
% A = struct('BF_ACC',-peak_timing_A_BF_m_A,'BF_TH',peak_timing_BF_TH_m_A,...
%     'TH_ACC',-peak_timing_A_TH_m_A,'TH_BF',-peak_timing_BF_TH_m_A);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_m_A;peak_timing_BF_TH_m_A;...
    -peak_timing_A_TH_m_A;-peak_timing_BF_TH_m_A};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('SLO, Morning');

subplot(324);
% A = struct('BF_ACC',-peak_timing_A_BF_n_A,'BF_TH',peak_timing_BF_TH_n_A,...
%     'TH_ACC',-peak_timing_A_TH_n_A,'TH_BF',-peak_timing_BF_TH_n_A);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_n_A;peak_timing_BF_TH_n_A;...
    -peak_timing_A_TH_n_A;-peak_timing_BF_TH_n_A};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('SLO, Evening');

subplot(325);
% A = struct('BF_ACC',-peak_timing_A_BF_m_B,'BF_TH',peak_timing_BF_TH_m_B,...
%     'TH_ACC',-peak_timing_A_TH_m_B,'TH_BF',-peak_timing_BF_TH_m_B);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_m_B;peak_timing_BF_TH_m_B;...
    -peak_timing_A_TH_m_B;-peak_timing_BF_TH_m_B};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('SLW, Morning');

subplot(326);
% A = struct('BF_ACC',-peak_timing_A_BF_n_B,'BF_TH',peak_timing_BF_TH_n_B,...
%     'TH_ACC',-peak_timing_A_TH_n_B,'TH_BF',-peak_timing_BF_TH_n_B);
% violin_plot(A,{'BF-ACC','BF-TH','TH-ACC','TH-BF'});
inputs = {-peak_timing_A_BF_n_B;peak_timing_BF_TH_n_B;...
    -peak_timing_A_TH_n_B;-peak_timing_BF_TH_n_B};
label = {'BF-ACC','BF-TH','TH-ACC','TH-BF'};
bar_plot_scatter(inputs,RGB,label);
ylabel({'Time peak';'(time frame)'});
title('SLW, Evening');

% healthy
figure(2);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,500,525]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(221);
plotbar_withstar({-peak_timing_A_BF_m_C,peak_timing_BF_TH_m_C},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title('morning , healthy control');

subplot(222);
plotbar_withstar({ - peak_timing_A_BF_n_C, peak_timing_BF_TH_n_C},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title(' night, healthy control');

subplot(223);
plotbar_withstar({-peak_timing_A_TH_m_C,-peak_timing_BF_TH_m_C},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('morning, healthy control');

subplot(224);
plotbar_withstar({ - peak_timing_A_TH_n_C,- peak_timing_BF_TH_n_C},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('night, healthy control');

% SLO
figure(3);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,500,525]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(221);
plotbar_withstar({-peak_timing_A_BF_m_A,peak_timing_BF_TH_m_A},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title('morning , SLO');

subplot(222);
plotbar_withstar({ - peak_timing_A_BF_n_A, peak_timing_BF_TH_n_A},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title(' night, SLO');

subplot(223);
plotbar_withstar({-peak_timing_A_TH_m_A,-peak_timing_BF_TH_m_A},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('morning, SLO');

subplot(224);
plotbar_withstar({ - peak_timing_A_TH_n_A,- peak_timing_BF_TH_n_A},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('night, SLO');

% SLW
figure(4);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,500,525]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(221);
plotbar_withstar({-peak_timing_A_BF_m_B,peak_timing_BF_TH_m_B},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title('morning , SLW');

subplot(222);
plotbar_withstar({ - peak_timing_A_BF_n_B, peak_timing_BF_TH_n_B},1000);
ylabel('peak timing');
xticklabels(['BF-ACC';'BF-TH ']);
title(' night, SLW');

subplot(223);
plotbar_withstar({-peak_timing_A_TH_m_B,-peak_timing_BF_TH_m_B},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('morning, SLW');

subplot(224);
plotbar_withstar({ - peak_timing_A_TH_n_B,- peak_timing_BF_TH_n_B},1000);
ylabel('peak timing');
xticklabels(['TH-ACC';'TH-BF ']);
title('night, SLW');

