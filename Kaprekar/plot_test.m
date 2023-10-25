% Pierce Zhang, CMOR220, Fall 2023, Project 4: Kaprekar's Constant
% Kaprekar.m
% Plot Kaprekar's quest for all valid four-digit integer starting values
% Last modified: 30 September 2023

% Driver
function plot_test
    data = quest(5,24);
    unique(data(:,24))
end

function [d] = dis(x)
    d = 0;
    for i = 0:x-1
        d = d + 10^i;
    end
end

function [n_next] = kaprekar_iter(n, x)
    string_num = num2str(x);
    if length(string_num) < n
        string_num = strcat('0',string_num);
    end
    asc = str2double(sort(string_num,'ascend'));
    desc = str2double(sort(string_num,'descend'));
    n_next = desc - asc;
end

% n is number of digits
function [data] = quest(n, maxiter)
    V = 10^(n-1):1:10^n-1;
    D = dis(n);
    V(mod(V,D) == 0) = [];

    dp = NaN(1,10^n);
    data = zeros(length(V),maxiter+1);
    data(:,1) = transpose(V);

    for i=1:maxiter
        for j=1:length(V)
            if (data(j,i) == 0)
                data(j,i+1) = 0;
            elseif (isnan(dp(data(j,i))))
                data(j,i+1) = kaprekar_iter(n, data(j,i));
                dp(data(j,i)) = data(j,i+1);
            else
                data(j,i+1) = dp(data(j,i));
            end
        end
    end
end