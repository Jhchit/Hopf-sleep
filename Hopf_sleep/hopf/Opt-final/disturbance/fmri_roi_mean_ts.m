function TS = fmri_roi_mean_ts(fmri,atlas, node_idx)
n_rois = length(node_idx);
[X,Y,Z,T] = size(fmri);
TS = zeros(T, n_rois);

% 设置滤波参数
fs = 1/3;
[b, a] = butter(4, [2*0.01/fs 2*0.1/fs], 'bandpass');

for roi = 1:n_rois
    roi_idx = node_idx(roi);

    mask = atlas==roi_idx;
    roi_dim = sum(mask, 'all');
    mask_coor = zeros(roi_dim, 3);

    if roi_dim > 0
        count = 1;
        for i = 1:X
            for j = 1:Y
                for k = 1:Z
                    if mask(i, j, k)
                        mask_coor(count, :) = [i j k];
                        count = count + 1;
                    end
                end
            end
        end

        ts = zeros(T, roi_dim);
        for i = 1:roi_dim
            idx = mask_coor(i, :);
            ts(:, i) = squeeze(fmri(idx(1), idx(2), idx(3), :));
        end
        ts = ts - mean(ts, 1);
        % 进行滤波操作
        ts = filtfilt(b, a, ts);
    else
        ts = 0;
    end
    TS(:, roi) = mean(ts, 2);
end
