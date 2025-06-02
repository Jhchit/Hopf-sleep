%%%% code for Fig. 3 and Fig.4
clear all;
close all;

%% data loading
filename = 'Hopf-Sleep\analysis\data for analysis'; % the name of the file containing empirical ALFF and FC
load([filename,'/alff_n.mat']);
load([filename,'/alff_m.mat']);
index_A = [1:18,55:59,65]; %SLO
index_B = [19:38,60:64,66,71:80];%SLW
index_C = [39:51,67:70];%HC
fc23 = load([filename,'/fc.mat']);
G_A = xlsread([filename,'/gender_age.xlsx']);
% regress out age and gender
    ACC_m = alff_m(1,:)'; 
    ACC_n = alff_n(1,:)'; 
    group = zeros(length(ACC_m),1);
    group(index_A) = 1;
    group(index_B) = 2;
    group(index_C) = 3;
    AGE = G_A(:,3);
    GENDER = G_A(:,2);
    tbl = table(AGE,GENDER,ACC_m,group);
    lme = fitlme(tbl,'ACC_m ~ AGE + GENDER + (1|group)');
    b1 = lme.Coefficients{2,2};
    b2 = lme.Coefficients{3,2};
    b0 = lme.Coefficients{1,2};
    ACC_m = ACC_m - (b1*AGE+b2*GENDER);
    
    tbl = table(AGE,GENDER,ACC_n,group);
    lme = fitlme(tbl,'ACC_n ~ AGE + GENDER + (1|group)');
    b1 = lme.Coefficients{2,2};
    b2 = lme.Coefficients{3,2};
    b0 = lme.Coefficients{1,2};
    ACC_n = ACC_n - (b1*AGE+b2*GENDER);

    gender_A = GENDER(index_A);
    gender_B = GENDER(index_B);
    gender_C = GENDER(index_C);
    age_A = AGE(index_A);
    age_B = AGE(index_B);
    age_C = AGE(index_C);
    fc23_m_A = fc23.A.m;
    fc23_n_A = fc23.A.n;
    fc23_m_B = fc23.B.m;
    fc23_n_B = fc23.B.n;
    fc23_m_C = fc23.C.m;
    fc23_n_C = fc23.C.n;
    fc23.A.m = regressout_age_gender(fc23.A.m',age_A,gender_A);
    fc23.A.n = regressout_age_gender(fc23.A.n',age_A,gender_A);
    fc23.B.m = regressout_age_gender(fc23.B.m',age_B,gender_B);
    fc23.B.n = regressout_age_gender(fc23.B.n',age_B,gender_B);
    fc23.C.m = regressout_age_gender(fc23.C.m',age_C,gender_C);
    fc23.C.n = regressout_age_gender(fc23.C.n',age_C,gender_C);





%  violin plot  for activation of ACC and FC TH-BF
flag = 1; % 1--- exclusing outliers， 2---doing nothing
var_num = 3;
figure(1);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(121);
A = struct('MH',fc23.A.m,'MSO',fc23.B.m,'MSW',fc23.C.m);
violin_plot(A,{'Healthy, morning','SLO, morning','SLW, morning'});
ylabel('FC');

subplot(122);
A = struct('MH',fc23.A.n,'MSO',fc23.B.n,'MSW',fc23.C.n);
violin_plot(A,{'Healthy, night','SLO, night','SLW, night'});
ylabel('FC');

figure(2);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(121);
plotbar_withstar({fc23.A.m,fc23.B.m,fc23.C.m},1000);
ylabel('FC');
xticklabels({'Healthy, morning','SLO, morning    ','SLW, morning    '});

subplot(122);
plotbar_withstar({fc23.A.n,fc23.B.n,fc23.C.n},1000);
ylabel('FC');
xticklabels({'Healthy, night','SLO, night    ','SLW, night    '});

figure(3);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(121);
if flag == 1
    alff_m_C = exclude_var(alff_m(1,index_C),var_num);
    alff_m_B = exclude_var(alff_m(1,index_B),var_num);
    alff_m_A = exclude_var(alff_m(1,index_A),var_num);
else
    alff_m_C = ACC_m(index_C);
    alff_m_B = ACC_m(index_B);
    alff_m_A = ACC_m(index_A);
end
A = struct('MH',alff_m_C,'MSO',alff_m_A,'MSW',alff_m_B);
violin_plot(A,{'Healthy, Morning','SLO, Morning','SLW, Morning'});
ylabel('ALFF');
title('Activation of ACC');

subplot(122);
if flag == 1
    alff_n_C = exclude_var(alff_n(1,index_C),var_num);
    alff_n_B = exclude_var(alff_n(1,index_B),var_num);
    alff_n_A = exclude_var(alff_n(1,index_A),var_num);
else
    alff_n_C = ACC_n(index_C);
    alff_n_B = ACC_n(index_B);
    alff_n_A = ACC_n(index_A);
end
A = struct('MH',alff_n_C,'MSO',alff_n_A,'MSW',alff_n_B);
violin_plot(A,{'Healthy, Night','SLO, Night','SLW, Night'});
ylabel('ALFF');
title('Activation of ACC');

figure(4);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度

subplot(121);
if flag == 1
    alff_m_C = exclude_var(alff_m(1,index_C),var_num);
    alff_m_B = exclude_var(alff_m(1,index_B),var_num);
    alff_m_A = exclude_var(alff_m(1,index_A),var_num);
else
    alff_m_C = ACC_m(index_C);
    alff_m_B = ACC_m(index_B);
    alff_m_A = ACC_m(index_A);
end
plotbar_withstar({alff_m_C,alff_m_A,alff_m_B},5000);
xticklabels(['Healthy, Morning';'SLO, Morning    ';'SLW, Morning    ']);
ylabel('ALFF');
title('Activation of ACC');

subplot(122);
if flag == 1
    alff_n_C = exclude_var(alff_n(1,index_C),var_num);
    alff_n_B = exclude_var(alff_n(1,index_B),var_num);
    alff_n_A = exclude_var(alff_n(1,index_A),var_num);
else
    alff_n_C = ACC_n(index_C);
    alff_n_B = ACC_n(index_B);
    alff_n_A = ACC_n(index_A);
end
plotbar_withstar({alff_n_C,alff_n_A,alff_n_B},5000)
xticklabels(['Healthy, Night';'SLO, Night    ';'SLW, Night    ']);
ylabel('ALFF');
title('Activation of ACC');

figure(5);
    set(gca, 'FontSize',8);
    set(gca,'FontName','Arial');
    set(gcf,'Position',[100,50,700,450]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
    set(gca,'Position',[.20 .17 .74 .60]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
    
    subplot(231); % HC
    ACC = alff_m(1,index_C)';
    BAS = alff_m(2,index_C)';
    THA = alff_m(3,index_C)';
    AGE = G_A(index_C,3);
    GENDER = G_A(index_C,2);
    FC23 = fc23.C.m;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('HC, Morning');
    ylabel('Regress Coef');
     %%%% stepwise regression- confirm
    fprintf('HC, Morning\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);

    subplot(234); % HC
    ACC = alff_n(1,index_C)';
    BAS = alff_n(2,index_C)';
    THA = alff_n(3,index_C)';
    AGE = G_A(index_C,3);
    GENDER = G_A(index_C,2);
    FC23 = fc23.C.n;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('HC, Evening');
    ylabel('Regress Coef');
    %%%% stepwise regression- confirm
    fprintf('HC,Evening\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);
   

    subplot(232); % SLO
    ACC = alff_m(1,index_A)';
    BAS = alff_m(2,index_A)';
    THA = alff_m(3,index_A)';
    AGE = G_A(index_A,3);
    GENDER = G_A(index_A,2);
    FC23 = fc23.A.m;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('SLO, Morning');
    ylabel('Regress Coef');
    %%%% stepwise regression- confirm
    fprintf('SLO, Morning\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);

    subplot(235); %SLO
    ACC = alff_n(1,index_A)';
    BAS = alff_n(2,index_A)';
    THA = alff_n(3,index_A)';
    AGE = G_A(index_A,3);
    GENDER = G_A(index_A,2);
    FC23 = fc23.A.n;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('SLO, Evening');
    ylabel('Regress Coef');
    %%%% stepwise regression- confirm
    fprintf('SLO, Evening\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);

    subplot(233); % SLW
    ACC = alff_m(1,index_B)';
    BAS = alff_m(2,index_B)';
    THA = alff_m(3,index_B)';
    AGE = G_A(index_B,3);
    GENDER = G_A(index_B,2);
    FC23 = fc23.B.m;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('SLW, Morning');
    ylabel('Regress Coef');
    %%%% stepwise regression- confirm
    fprintf('SLW, Morning\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);

    subplot(236); % SLW
    ACC = alff_n(1,index_B)';
    BAS = alff_n(2,index_B)';
    THA = alff_n(3,index_B)';
    AGE = G_A(index_B,3);
    GENDER = G_A(index_B,2);
    FC23 = fc23.B.n;
    tbl = table(ACC,AGE,GENDER,FC23);
    lme = fitlme(tbl,'ACC ~ AGE+GENDER+FC23');
    plot_regression_bar(lme,[1,2,3],{'Age','Gender','FC'});
    title('SLW, Evening');
    ylabel('Regress Coef');
    %%%% stepwise regression- confirm
    fprintf('SLW, Evening\n',repmat('-', 1, 35));
    stlm = stepwiselm(tbl,'ACC ~ AGE+GENDER+FC23','Upper','linear','Verbose', 0);
    disp(stlm.Coefficients);


%% functions

function [] = plot_regression_bar(lme,idx,names)
 % lme --- the model regressed
 % idx --- the index for the regressors
 % names --- the names of presented regressors
     if length(names) ~= length(idx)
         error('wrong input: names and idx do not match');
     end

        y = table2array(dataset2table(lme.Coefficients(:,2))); 
        y= y(1+idx);

        neg = table2array(dataset2table(lme.Coefficients(:,7))); 
        neg = neg(1+idx);
        
        pos = table2array(dataset2table(lme.Coefficients(:,8)));  
        pos = pos(1+idx);

        RGB1 = [0.80392,0.78824,0.78824]; % light grey[0.28235,0.23922,0.5451]; % highly infiltrated regions[0.86667,0.62745,0.86667]; % light purple[0.80392,0.78824,0.78824]; % light grey[0.5451,0.53725,0.53725]; % dark grey[0.80392,0.78824,0.78824]; % light grey for low-grade%[0.00, 0.45, 0.74];%[0.00,0.45,0.74];%蓝色  [217,205,204]/264;
        RGB2 = [0.5451,0.53725,0.53725]; % dark grey[0.8549,0.64706,0.12549]; % lowly infiltrated regions[0.84706,0.74902,0.84706]; % light red[1,0.89412,0.7098];% healthy%[0.85, 0.33, 0.10];%[0.85,0.33,0.10];%红色  [122,150,171]/264;
        RGB3 = [0.91373,0.58824,0.47843]; % surrounding regions[1,0.71373,0.75686]; % light pink [0.91373,0.58824,0.47843]; % surrounding regions[0.5451,0.53725,0.53725]; % dark grey high-grade [0.82353,0.41176,0.11765]; % surrounding regions[191,141,122]/264;%[0.9290 0.6940 0.1250];
        
        RGB = [RGB1;RGB2;RGB3];
        
        m = size(y,1);
        n = size(y,2);
        x = 1 : m;
        h = bar(x, y,'BarWidth',0.8);
        for i = 1 : m
            for j = 1 : n
                h(1, j).FaceColor = 'flat';
                h(1, j).CData(i,:) = RGB(i,:);
                h(1,j).EdgeColor = 'flat';
            end
        end
        % 绘制误差线
        hold on
        errorbar(x, y, y-neg, pos-y, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
        hold off
        
        names = names';
        set(gca, 'XTickLabel', names,'FontSize',12,'FontWeight','bold','XTickLabelRotation',-45);


        pvalue1 = table2array(dataset2table(lme.Coefficients(:,6)));  
        pvalue = pvalue1(idx+1);
        idx1= find(pvalue < 0.05);
end

function [] = violin_plot(inputs,label)
% inputs --- struct, M*1

vector2plot = inputs;
B = fieldnames(inputs);
violins = violinplot(vector2plot);
colors = cbrewer2('seq','YlGnBu',size(B,1),'linear');
for idx = 1:size(B,1)
    violins(idx).ViolinAlpha = 0.7;
    violins(idx).EdgeColor = [0 0 0];
    violins(idx).BoxColor = [0 0 0];
    violins(idx).ViolinColor = colors(idx,:);
    violins(idx).ScatterPlot.MarkerFaceColor = [0 0 0]; 
    violins(idx).ScatterPlot.MarkerFaceAlpha = 1; 
end
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 2, 'TickDir', 'out'); 
title('', 'FontSize', 12, 'FontName', 'Arial');
xticklabels(label);

end

function data_screened = exclude_var(data,var_num)
% data --- vector
% var_num --- how many times of varance should be included into the output
% data/-screened
  std_data = std(data)/sqrt(length(data));
  avg_data = mean(data);
  data_screened = [];
  num = 0;
  for i = 1:length(data)
      if data(i) <= avg_data + var_num*std_data & data(i) >= avg_data - var_num*std_data
         num = num + 1;
         data_screened(num,1) = data(i);
      end
  end
end

function data_regress = regressout_age_gender(data,age,gender)
    tbl = table(age,gender,data);
    lme = fitlme(tbl,'data ~ age + gender');
    b1 = lme.Coefficients{2,2};
    b2 = lme.Coefficients{3,2};
    b0 = lme.Coefficients{1,2};
    data_regress = data- (b1*age+b2*gender);
end
