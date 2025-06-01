%左右脑区取平均计算
clc;
clear;

load('corridatotal.mat','out','Cfg','SC_cell','Hopf');%载入hopf模型的拟合结果

OUT = out;
Rta = out.Rta;
Rta(91) = 0.001;   %31是acc左脑区%91是下丘脑右
Rta(92) = 0.001;   %32是acc右脑区%92是下丘脑左
n = 19;

result(100) = struct('Alff',[],'mean_alff',[],'FSim_x',[],'meanfc_x',[]);
meanfc_y = zeros(4,100);
FSim_y = zeros(4,100);
for i = 1:100
    out = OUT(i);%选取指定out值

    %Alff是是熬夜睡眠无障碍组的acc的alff， mean_alff是平均alff
    [Alff,mean_alff,meanfc,FSim,xs] = a_change(out,Cfg,SC_cell,Hopf,Rta,n);%这里的meanfc指的是一组人的平均，不是100次扰动的平均%FSim是模拟的扰动前的FC
    
    %干扰后
    xs(93,:) = mean(xs([21,41],:),1);%基底前脑左
    xs(94,:) = mean(xs([22,42],:),1);%基底前脑右

    %干扰前
    XS = out.xs;
    XS(93,:) = mean(XS([21,41],:),1);
    XS(94,:) = mean(XS([22,42],:),1);

    xs(95,:) = mean(xs(31:32,:),1);%acc
    xs(96,:) = mean(xs(93:94,:),1);%基底前脑
    xs(97,:) = mean(xs(91:92,:),1);%下丘脑
    meanfc = observable_FC(xs);
    
    XS(95,:) = mean(XS(31:32,:),1);%acc
    XS(96,:) = mean(XS(93:94,:),1);%基底前脑
    XS(97,:) = mean(XS(91:92,:),1);%下丘脑
    FSIM = observable_FC(XS);


    meanfc_x = meanfc(95:97,95:97);%提取扰动后acc与下丘脑的相关矩阵
    FSim_x = FSIM(95:97,95:97);   %扰动前acc与下丘脑的相关矩阵



    result(i).Alff = Alff;
    result(i).mean_alff = mean_alff;
    result(i).FSim_x = FSim_x;
    result(i).meanfc_x = meanfc_x;    
end


