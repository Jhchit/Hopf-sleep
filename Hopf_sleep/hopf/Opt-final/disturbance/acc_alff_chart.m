arr_fc = zeros(1,100);%   Rta(31) = 4;
for i = 1:100
    arr_fc(i) = result(i).meanfc_x(2,2);
end

% 排序并获取排序后的索引
[sortedArr, originalIndices] = sort(arr_fc);

% 显示排序后的数组及原始索引
disp('排序后的数组:');
disp(sortedArr);

disp('排序前数的序号（按排序后次序排列）:');
disp(originalIndices);

arr_alff = zeros(1,100);
for i = 1:100
    arr_alff(i) = result(originalIndices(i)).mean_alff{2}(2);
end

%散点图加线性回归
scatter(sortedArr,arr_alff,5,'r','filled');
p = polyfit(sortedArr,arr_alff,1);
y_fit = polyval(p,sortedArr);
hold on
plot(sortedArr,y_fit);
xlabel('acc(右)与下丘脑(左)的FC');
ylabel('acc(右)_alff');



%扰动前后，acc的alff的变化
accsim_alff = zeros(2,100);
accr_alff = zeros(2,100);
for i = 1:100
    accsim_alff(1,i) = result(i).mean_alff{1}(1);
    accsim_alff(2,i) = result(i).mean_alff{1}(2);
    accr_alff(1,i) = result(i).mean_alff{2}(1);
    accr_alff(2,i) = result(i).mean_alff{2}(2);
end

%箱型图
for i = 1:2
    subplot(1,2,i)
    group = [repmat('干扰前', size(accsim_alff(i,:),2), 1); repmat('干扰后', size(accr_alff(i,:),2), 1)];
    boxplot([accsim_alff(i,:)' accr_alff(i,:)'],group);
     ylabel('组平均acc的激活程度(alff)');
     if i == 1
        title('acc(左)');
    end
    if i == 2
        title('acc(右)');
    end
end

%平均水平，柱状图
meandata_1 = mean(accsim_alff,2);
meandata_2 = mean(accr_alff,2);
subplot(1,2,1);
bar([meandata_1(1) meandata_2(1)]);
subplot(1,2,2);
bar([meandata_1(2) meandata_2(2)]);


