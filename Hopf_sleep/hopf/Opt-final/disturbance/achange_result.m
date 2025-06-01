%左右脑分开分析
clc;
clear;

load('corridatotal.mat','out','Cfg','SC_cell','Hopf');%载入hopf模型的拟合结果

OUT = out;
Rta = out.Rta;

%改变分叉参数
Rta(91) = -0.01;   %31是acc左脑区%91是下丘脑右
Rta(92) = -0.01;  %32是acc右脑区%92是下丘脑左


n = 19;%早上填18，晚上19

result(100) = struct('Alff',[],'mean_alff',[],'FSim_x',[],'meanfc_x',[]);
meanfc_y = zeros(4,100);
FSim_y = zeros(4,100);
for i = 1:100
    out = OUT(i);%选取指定out值
    [Alff,mean_alff,meanfc,FSim,xs] = a_change(out,Cfg,SC_cell,Hopf,Rta,n);%这里的meanfc指的是一组人的平均，不是100次扰动的平均%FSim是模拟的扰动前的FC
    
    %干扰后
    xs(93,:) = mean(xs([21,41],:),1);
    xs(94,:) = mean(xs([22,42],:),1);

    %干扰前
    XS = out.xs;
    XS(93,:) = mean(XS([21,41],:),1);
    XS(94,:) = mean(XS([22,42],:),1);
    FSIM = observable_FC(XS);

    meanfc_x = meanfc([31:32,91:94],[31:32,91:94]);%提取acc与下丘脑的相关矩阵
    FSim_x = FSIM([31:32,91:94],[31:32,91:94]);   


    result(i).Alff = Alff;
    result(i).mean_alff = mean_alff;
    result(i).FSim_x = FSim_x;
    result(i).meanfc_x = meanfc_x;    
end


%散点图
% i =1;
% subplot(2,2,i);
% scatter(1:100,FSim_y(i,:),5,'r','filled');
% hold on
% scatter(1:100,meanfc_y(i,:),5,'b','filled');
% legend('不干扰','干扰使acc激活程度下降');
% ylabel('组平均fc');
% xlabel('模拟次数');
% title('熬夜不失眠组,干扰acc激活程度前后对比,后acc(左)与下丘脑(右)');
% 
% i =2;
% subplot(2,2,i);
% scatter(1:100,FSim_y(i,:),5,'r','filled');
% hold on
% scatter(1:100,meanfc_y(i,:),5,'b','filled');
% legend('不干扰','干扰使acc激活程度下降');
% ylabel('组平均fc');
% xlabel('模拟次数');
% title('熬夜不失眠组,干扰acc激活程度前后对比,acc(左)与下丘脑(左)');
% 
% i =3;
% subplot(2,2,i);
% scatter(1:100,FSim_y(i,:),5,'r','filled');
% hold on
% scatter(1:100,meanfc_y(i,:),5,'b','filled');
% legend('不干扰','干扰使acc激活程度下降');
% ylabel('组平均fc');
% xlabel('模拟次数');
% title('熬夜不失眠组,干扰acc激活程度前后对比,acc(右)与下丘脑(右)');
% 
% i =4;
% subplot(2,2,i);
% scatter(1:100,FSim_y(i,:),5,'r','filled');
% hold on
% scatter(1:100,meanfc_y(i,:),5,'b','filled');
% legend('不干扰','干扰使acc激活程度下降');
% ylabel('组平均fc');
% xlabel('模拟次数');
% title('熬夜不失眠组,干扰acc激活程度前后对比,acc(右)与下丘脑(左)');


% %箱型图
% for i= 1:4
% 
%     subplot(2,2,i);
%     group = [repmat('干扰前', size(FSim_y(i,:),2), 1); repmat('干扰后', size(meanfc_y(i,:),2), 1)];
%     boxplot([FSim_y(i,:)' meanfc_y(i,:)'],group);
%     ylabel('组平均fc');
%    
%     if i == 1
%         title('熬夜不失眠组,干扰下丘脑激活程度前后对比,acc(左)与下丘脑(右)');
%     end
%     if i == 2
%         title('熬夜不失眠组,干扰下丘脑激活程度前后对比,acc(左)与下丘脑(左)');
%     end
%     if i == 3
%         title('熬夜不失眠组,干扰下丘脑激活程度前后对比,acc(右)与下丘脑(右)');
%     end
%     if i == 4
%         title('熬夜不失眠组,干扰下丘脑激活程度前后对比,acc(右)与下丘脑(左)');
%     end    
% end



% %绘制平均水平柱状图。并增加显著标志
% i = 1;
% p1 = 0.002;
% data_mean = [mean(FSim_y(i,:)),mean(meanfc_y(i,:))];   %此时m=1,n=2
% data_std = [std(FSim_y(i,:)),std(meanfc_y(i,:))];
% 
% 
% % data_mean = [mean(accsim_alff(i,:)),mean(accr_alff(i,:))]; 该数据计算acc的alff
% % data_std = [std(accsim_alff(i,:)),std(accr_alff(i,:))];
% 
% % 绘制条形图
% figure;
% y = data_mean;
% neg = data_std;
% pos = data_std;
% m = 1;
% n = size(y,2);
% name = categorical({'干扰下丘脑激活程度前后对比,acc(右)与下丘脑(左)'});
% h = bar(name,y,0.4);
% title('熬夜不失眠组');
% legend('干扰前','干扰后');
% ylabel('100次扰动下的平均fc')
% 
% 
% 
% % 获取误差线 x 值
% xx = zeros(m, n);
% for i = 1 : n
%     xx(:, i) = h(1, i).XEndPoints';
% end
% 
% % 绘制误差线
% hold on
% errorbar(xx, y, neg, pos, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
% 
% %绘制显著性差异
% t1 = 1; %第1组
% n1 = 1; %第1组的第n1个柱子
% n2 = 2; %第1组的第n2个柱子
% %获取显著性差异线的x坐标
% x1 = xx(t1,n1); x2 = xx(t1,n2);
% %获取显著性差异线的y坐标
% ySig = max(y(t1,n1)+pos(t1,n1),y(t1,n2)+pos(t1,n2));
% %绘制显著性差异
% sigline([x1,x2],ySig, p1); % p1是上面T检验得到的p值

% sum1 = 0;
% sum2 = 0;
% for i = 1:100
%     sum1 = sum1 + result(i).mean_alff{1};
%     sum2 = sum2 + result(i).mean_alff{2};
% end
% ans1 = sum1/100;
% ans2 = sum2/100;
% 
% 
% for i = 1:100
%     n1(i) =result(i).mean_alff{1};
%     n2(i) =result(i).mean_alff{2};
% end
% 
% 
% 
% 
% fc1 = zeros(1,100);
% for i = 1:100
%     fc1(i) = mean(mean(result(i).FSim_x(3:4,5:6)));
% end
% meanfc1 = mean(fc1);
% fc2 = zeros(1,100);
% for i = 1:100
%     fc2(i) = mean(mean(result(i).meanfc_x(3:4,5:6)));
% end
% 
% meanfc2 = mean(fc2);

