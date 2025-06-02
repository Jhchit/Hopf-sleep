function [r_lag] = Pearson_time_lage(sig1,sig2,max_lag)
%  sig1 --- N*1, sig2 --- N*1, max_lag --- the largest time lag considered
%  output: r_lag --- 1*N, -max_lag ~ 0 ~ +max_lag
N = length(sig1);

 
 for lag = 0:max_lag
     sig1_1 = sig1(1+lag:N);
     sig2_1 = sig2(1:length(sig1_1));
     r = corrcoef(sig1_1,sig2_1);
     [L,W] = size(r);
     if L > 1
        r_lag1(1,lag+1) = r(1,2);
     else if L == 1
         r_lag1(1,lag+1) = r;
     else 
         break
     end
     end

 end
 for lag = 1:max_lag
     sig2_2 = sig2(1+lag:N);
     sig1_2 = sig1(1:length(sig2_2));
     r = corrcoef(sig1_2,sig2_2);
     [L,W] = size(r);
     if L > 1
        r_lag2(1,max_lag-lag+1) = r(1,2);
     else if L == 1
             r_lag2(1,max_lag-lag+1) = r;
     else 
         break
     end
     end     
 end
r_lag = [r_lag2,r_lag1];
end