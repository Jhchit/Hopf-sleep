%meanfc:干扰后的平均fc矩阵;xs是干扰后的时序；FSim是out里的模拟FC
function [Alff,mean_alff,meanfc,FSim,xs] = a_change(out,Cfg,SC_cell,Hopf,Rta,n)
new_a = Rta;
new_G = out.RtG;
SC_xs = SC_cell{1};
TRsec = Cfg.filt.TRsec;
nNodes = Cfg.nNodes;
long_Total = Hopf.long_Total(1);
w = Hopf.w{1};
val = Hopf.val{1};
Gmethod = Hopf.Gmethod;
xs = resim_Hopf(new_a,new_G,SC_xs,TRsec,nNodes,long_Total,w,val,Gmethod);

if Cfg.filt.bpass==1 

    xs = filtroign(xs,Cfg.filt.TRsec,Cfg.filt.lb,Cfg.filt.ub); 

end


T = 60; % 时间点数目%时间点数随n而变
N = 1;%roi数量;%左右半脑的acc若取平均，N为1；反之为2


fs = 1 / 3; % TR是重复时间，也就是两个连续采样点之间的时间间隔


%提取出acc的时序
xs_acc = cell(1,2);
if N == 1
    xs_acc{1} = mean(out.xs(31:32,:));%干扰前acc的时序
    xs_acc{2} = mean(xs(31:32,:));%干扰后acc的时序
else
    xs_acc{1} = out.xs(31:32,:);
    xs_acc{2} = xs(31:32,:);
end

Alff = cell(1,2);
for x = 1:2
    xs_accN = xs_acc{x};

    %19人是否分开算
    if n == 19
        ts_1 = cell(1,n);
        for i =1:n
            ts_1{i} = xs_accN(:,T*(i-1)+1:T*i);
        end
    else
        ts_1 = {xs_accN};
    end
   
    alff = zeros(N, n);
    
    %计算alff
    for j = 1:size(alff,2)
        ts = (ts_1{j});
        
        for i = 1:N
            % 傅里叶变换
            fft_ts = fft(ts(i, :)); 
            
            % 计算频谱
            P2 = abs(fft_ts / T); % 双边谱
            P1 = P2(1:T/2+1);
            P1(2:end-1) = 2 * P1(2:end-1); % 单边谱
            
            % 计算ALFF
            low_freq_range = 0.01; % 低频范围（例如0.01-0.1 Hz）
            high_freq_range = 0.1;
            low_freq_idx = round(low_freq_range * T * fs) + 1; % +1 是因为 MATLAB 索引从 1 开始
            high_freq_idx = round(high_freq_range * T * fs) + 1;
            
            alff(i,j) = sum(P1(low_freq_idx:high_freq_idx));
        end
    end
    Alff{x} = alff;

mean_alff{x} = mean(alff,2);

end



meanfc = observable_FC(xs);%干扰后的平均fc矩阵

FSim = out.FSim;

% %计算每个人的fc矩阵
% fc1 = cell(1,n);%干扰前
% fc2 = cell(1,n);%干扰后
% for i =1:n
%     Ts = out.xs(:,T*(i-1)+1:T*i);
%     fc = corr(Ts');
%     fc1{i} = fc(31:32,91:92);
%     Ts = xs(:,T*(i-1)+1:T*i);
%     fc = corr(Ts');
%     fc2{i} = fc(31:32,91:92);
% end


% 
% for i = 1:n
%     fc1_1(i) = fc1{i}(1);
%     fc2_1(i) = fc2{i}(1);
% end

end