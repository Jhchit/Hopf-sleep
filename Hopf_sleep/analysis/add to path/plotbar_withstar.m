function [ pvalue ] = plotbar_withstar(input,permutation_num)
        % input--- 1*M cell, each component contains N dots inside a
        % group/bar, M denotes groups/bars
        % p --- 1* (M*(M-1)/2)
        x = 1:length(input);
        pvalue = [];
        
        for i = 1:length(input)
            A = input{i};
            A(find(isnan(A) == 1)) = [];
            data(1,i) = mean(A);
            err(1,i)  = std(A);
            errhigh(1,i) = data(1,i)-(std(A)/sqrt(length(A)));
            errlow(1,i)  = data(1,i)+(std(A)/sqrt(length(A)));
           
        end
        bar(x,data);
        hold on
        er = errorbar(x,data,errlow-data);    
        er.Color = [0, 0, 0];                            
        er.LineStyle = 'none'; 
        yMargin = 0.5;  % 设置边距，可以根据实际需求调整
        yMin = min(errhigh) * (1+yMargin);
        yMax = max(errlow) * (1+yMargin);
        num = 0;
        for i = 1:(length(input)-1)
            if i ~= length(input)
                for j = (i+1):length(input)
                    num = num+1;
                    inp1 = input{i};
                    inp2 = input{j};
                   
                    [p, observeddifference, effectsiz] = permutationTest(inp1, ...
                                                    inp2, permutation_num);
%                     df = length(inp1)+length(inp2)-2;
%                     mu    = [mean(inp1); mean(inp2)];
%                     sigma = [std(inp1)/sqrt(length(inp1)); std(inp2)/sqrt(length(inp2))];
%                     %calculate t-statistic
%         
%                     stat = abs(diff(mu)) / sqrt(sigma(1)^2+sigma(2)^2);
%                     df = (sigma(1)^2+sigma(2)^2)^2/(sigma(1)^4/(length(inp1)-1) + sigma(2)^4/(length(inp2)-1));
%                     p = 2*tcdf(-abs(stat), df);

                    pvalue = [pvalue,p];
                    txt = num2str(p);
                    %{
                    if p > 0.05
                    txt = 'N.A';
                    else if p <= 0.05 & p > 0.001
                            txt = '*';
                    else if p <= 0.001 & p > 0.0001
                            txt = '**';
                    else if p <= 0.0001
                        txt = '***';
                    end
                    end
                    end
                    end
                    %}
                    
                    if mean(inp1) > 0 & mean(inp2) > 0
                        text(mean([i,j]),max(errlow)*(1.25+num*0.1),[txt],...
                            FontSize=10,FontWeight='bold');
                        %'difference:',num2str(mean(inp1)-mean(inp2))
                        line([i,j],[max(errlow)*(1.1+num*0.1),max(errlow)*(1.1+num*0.1)]);
%                         ylim(sort([0 max(errlow)*2]));
                    else if mean(inp1) < 0 & mean(inp2) < 0
                            height = min([mean(inp1),mean(inp2)]);
                            text(mean([i,j]),height*(1.25+num*0.1),[txt],....
                                FontSize=10,FontWeight='bold');
                            line([i,j],[height*(1.1+num*0.1),height*(1.1+num*0.1)]);
%                             ylim(sort([0 min(errhigh)*2]));
                    else
                        height = max([mean(inp1),mean(inp2)]);
                        text(mean([i,j]),height*(1.4+num*0.1),[txt],...
                            FontSize=10,FontWeight='bold');
                        line([i,j],[height*(1.2+num*0.1),height*(1.2+num*0.1)]);
%                         ylim(sort([min(errhigh)*2 max(errlow)*2]));
                    end
                    end
                end
            end
        end
        
        for i = 1:length(input)
            A = input{i};
            A(find(isnan(A) == 1)) = [];
            avg(1,i) = mean(A)
            text(i,avg(i)*1.2,num2str(avg(i)))
        end 
            
            
        hold off;
        ylim([yMin, yMax]);
        box off;
end

