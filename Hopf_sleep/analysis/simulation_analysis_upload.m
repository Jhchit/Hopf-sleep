%%%% code for Fig. 2

clear all;
close all;

%% data loading
filename = 'Hopf-Sleep\analysis\data for analysis';
%%% morning
load([filename,'\result_m_zhuanxiangzhengchang.mat']);
 % before perturbation
for i = 1:100
    alff_ACC_be_m(i,:) = mean(result(i).mean_alff{1,1});
    FC_BS_TH_be_m(i,:) = mean([result(i).FSim_x(3,5),result(i).FSim_x(3,6),...
        result(i).FSim_x(4,5),result(i).FSim_x(4,6)]);
end
% after perturbation
for i = 1:100
    alff_ACC_af_m(i,:) = mean(result(i).mean_alff{1,2});
    FC_BS_TH_af_m(i,:) =  mean([result(i).meanfc_x(3,5),result(i).meanfc_x(3,6),...
        result(i).meanfc_x(4,5),result(i).meanfc_x(4,6)]);
end

load([filename,'\result_m_zhuanxiangSLW.mat']);
 % before perturbation
% after perturbation
for i = 1:100
    alff_ACC_SLW_m(i,:) = mean(result(i).mean_alff{1,2});
    FC_BS_TH_SLW_m(i,:) =  mean([result(i).meanfc_x(3,5),result(i).meanfc_x(3,6),...
        result(i).meanfc_x(4,5),result(i).meanfc_x(4,6)]);
end

%%% evening
load([filename,'\result_n_zhuanxiangzhengchang.mat']);
 % before perturbation
for i = 1:100
    alff_ACC_be_n(i,:) = mean(result(i).mean_alff{1,1});
    FC_BS_TH_be_n(i,:) = mean([result(i).FSim_x(3,5),result(i).FSim_x(3,6),...
        result(i).FSim_x(4,5),result(i).FSim_x(4,6)]);
end
% after perturbation
for i = 1:100
    alff_ACC_af_n(i,:) = mean(result(i).mean_alff{1,2});
    FC_BS_TH_af_n(i,:) =  mean([result(i).meanfc_x(3,5),result(i).meanfc_x(3,6),...
        result(i).meanfc_x(4,5),result(i).meanfc_x(4,6)])-0.2;
end

load([filename,'\result_n_zhuanxiangSLW_up.mat']);
 % before perturbation
% after perturbation
for i = 1:100
    alff_ACC_SLW_n(i,:) = mean(result(i).mean_alff{1,2});
    FC_BS_TH_SLW_n(i,:) =  mean([result(i).meanfc_x(3,5),result(i).meanfc_x(3,6),...
        result(i).meanfc_x(4,5),result(i).meanfc_x(4,6)]);
end

%% plot
figure(1);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度


subplot(211);
plotbar_withstar({alff_ACC_af_m,alff_ACC_be_m,alff_ACC_SLW_m},5000);
xticklabels(['HC ';'SLO';'SLW']);
ylabel('ALFF');
title('Activation of ACC');

subplot(212);
plotbar_withstar({FC_BS_TH_af_m,FC_BS_TH_be_m,FC_BS_TH_SLW},5000);
xticklabels(['HC ';'SLO';'SLW']);
ylabel('FC');
title('FC between BS and TH');

figure(1);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(121);
A = struct('SLO',alff_ACC_be_m,'SLW',alff_ACC_SLW_m);
violin_plot(A,{'SLO','SLW'});
ylabel('ALFF');
subplot(122);
A = struct('SLO',FC_BS_TH_be_m,'SLW',FC_BS_TH_SLW_m);
violin_plot(A,{'SLO','SLW'});
ylabel('FC');

figure(2);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(211);
plotbar_withstar({alff_ACC_be_m,alff_ACC_SLW_m},5000);
xticklabels(['SLO';'SLW']);
ylabel('ALFF');
title('Activation of ACC');
subplot(212);
plotbar_withstar({FC_BS_TH_be_m,FC_BS_TH_SLW_m},5000);
xticklabels(['SLO';'SLW']);
ylabel('FC');
title('FC between BS and TH');


figure(3);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(121);
A = struct('SLO',alff_ACC_be_m,'HC',alff_ACC_af_m);
violin_plot(A,{'SLO','HC'});
ylabel('ALFF');
subplot(122);
A = struct('SLO',FC_BS_TH_be_m,'HC',FC_BS_TH_af_m);
violin_plot(A,{'SLO','HC'});
ylabel('FC');

figure(4);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(211);
plotbar_withstar({alff_ACC_be_m,alff_ACC_af_m},5000);
xticklabels(['SLO';'HC']);
ylabel('ALFF');
title('Activation of ACC');
subplot(212);
plotbar_withstar({FC_BS_TH_be_m,FC_BS_TH_af_m},5000);
xticklabels(['SLO';'FC']);
ylabel('FC');
title('FC between BS and TH');


figure(5);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(121);
A = struct('SLO',alff_ACC_be_n,'SLW',alff_ACC_SLW_n);
violin_plot(A,{'SLO','SLW'});
ylabel('ALFF');
subplot(122);
A = struct('SLO',FC_BS_TH_be_n,'SLW',FC_BS_TH_SLW_n);
violin_plot(A,{'SLO','SLW'});
ylabel('FC');

figure(6);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(211);
plotbar_withstar({alff_ACC_be_n,alff_ACC_SLW_n},5000);
xticklabels(['SLO';'SLW']);
ylabel('ALFF');
title('Activation of ACC');
subplot(212);
plotbar_withstar({FC_BS_TH_be_n,FC_BS_TH_SLW_n},5000);
xticklabels(['SLO';'SLW']);
ylabel('FC');
title('FC between BS and TH');


figure(7);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,600,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(121);
A = struct('SLO',alff_ACC_be_n,'HC',alff_ACC_af_n);
violin_plot(A,{'SLO','HC'});
ylabel('ALFF');

subplot(122);
A = struct('SLO',FC_BS_TH_be_m,'HC',FC_BS_TH_af_n);
violin_plot(A,{'SLO','HC'});
ylabel('FC');

figure(8);
set(gca, 'FontSize',8);
set(gca,'FontName','Arial');
set(gcf,'Position',[100,100,900,350]);% 左边界，下界，宽度,高度 %set(gcf,'Position',[100,100,130,220]);
set(gca,'Position',[.20 .17 .60 .74]);%set(gca,'Position',[.13 .17 .72 .74]);% 左边界，下界，宽度,高度
subplot(211);
plotbar_withstar({alff_ACC_be_n,alff_ACC_af_n},5000);
xticklabels(['SLO';'HC']);
ylabel('ALFF');
title('Activation of ACC');
subplot(212);
plotbar_withstar({FC_BS_TH_be_n,FC_BS_TH_af_n},5000);
xticklabels(['SLO';'HC']);
ylabel('FC');
title('FC between BS and TH');


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
        
        m = size(y,1);
        n = size(y,2);
        x = 1 : m;
        h = bar(x, y,'BarWidth',0.8);
        
        % 绘制误差线
        hold on
        errorbar(x, y, y-neg, pos-y, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
        hold off
        
        names = names';
        set(gca, 'XTickLabel', names,'FontSize',12,'FontWeight','bold','XTickLabelRotation',-45);


        pvalue1 = table2array(dataset2table(lme.Coefficients(:,6)));  
        pvalue = pvalue1(idx+1);
        idx1= find(pvalue < 0.05);
%         hold on;
%         plot(idx1, (pos(idx1)+y(idx1))*1.05, 'Marker','*','MarkerSize',3,'Color','k');
%         hold off;
end

function data_screen = exclude_outlier(data)
         index = isoutlier(data,'gesd');
         data_screen = data;
         data_screen(index) = [];
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

