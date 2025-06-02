function [] = bar_plot_scatter(inputs,RGB,label)
% input --- cell,M*1, M labels/groups
% RGB --- matrix, M*3, each row denotes one color
if size(inputs,1) ~= size(RGB,1)
    error('wrong input or RGB');
end
%% 数据 % m * n的矩阵，m对应的是一共有几组，n对应的是每组有几个柱子
for i = 1:length(inputs)
    A = sign(inputs{i,1}).*log(abs(inputs{i,1}));
    A(find(isnan(A)==1)) = 0;
    inputs{i,1} = A;
    data_mean(i,1) = mean(inputs{i,1});
    data_std(i,1) = std(inputs{i,1});
end


y = data_mean; 
% neg = data_std(:,[1,3]); 
% pos = data_std(:,[2,4]); 
neg = data_std; 
pos = data_std; 
m = size(y,1);
n = size(y,2);
x = 1 : m;
h = bar(x, y,'BarWidth',0.6);
%% 单独设置第i个组第j个柱子的颜色
for i = 1 : m
    for j = 1 : n
        h(1, j).FaceColor = 'flat';
        h(1, j).CData(i,:) = RGB(i,:);
        h(1,j).EdgeColor = 'flat';
    end
end
%% 获取误差线 x 值
xx = zeros(m, n);
for i = 1 : n
    xx(:, i) = h(1, i).XEndPoints';
end

% 绘制误差线
hold on
errorbar(xx, y, neg, pos, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
% errorbar(xx, y, y-neg, pos-y, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
hold off

% scatter
hold on;
for i = 1:length(inputs)
    x = i*ones(size(inputs{i,1}));
    scatter(x,inputs{i,1},10,'filled','b');
end
hold off;
yMin = min(data_mean-data_std);
yMax = max(data_mean+data_std);
ylim([yMin, yMax]*1.5);
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'LineWidth', 2, 'TickDir', 'out'); 
title('', 'FontSize', 12, 'FontName', 'Arial');
xticklabels(label);

end