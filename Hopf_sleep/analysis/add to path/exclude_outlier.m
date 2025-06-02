function data_screen = exclude_outlier(data,excluding_method)
% excluding_method --- struct,1*2; 1st elemnt --- use isoutlier or use
% standard deviation, 2nd element --- if is outlier, the specific method
% for this function, if standard devaition, how many times of STD are used.
         if excluding_method.method(1) == 'i' % is
            index = isoutlier(data,excluding_method.para);
         elseif excluding_method.method(1)  == 's' % std
             avg = mean(data);
             STAD = std(data);
             times = excluding_method.para;
             index = find((data<= avg - times*STAD) | (data>= avg + times*STAD));
         else 
             error('wrong input for excluding_method');
         end
         data_screen = data;
         data_screen(index) = [];
end