% Pierce Zhang, CMOR220, Fall 2023, Project 4: Kaprekar's Constant
% Kaprekar.m
% Plot Kaprekar's quest for all valid four-digit integer starting values
% Last modified: 30 September 2023

% Driver
function Kaprekar
    figure();
    plot(0:12,quest(5,12),'-o');
    title("Kaprekar Quest for 4-digit Values");
    ylabel("Values"); xlabel("Iteration Count");
    grid on;
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
    seen = zeros(1,10^n);
    data = zeros(length(V),maxiter+1);
    data(:,1) = transpose(V);

    for i=1:maxiter
        seen(:) = 0;
        for j=1:length(V)
            if isnan(data(j,i))
                data(j,i+1) = NaN;
            elseif seen(data(j,i)) == 1
                data(j,i+1) = NaN;
            elseif isnan(dp(data(j,i)))
                data(j,i+1) = kaprekar_iter(n, data(j,i));
                dp(data(j,i)) = data(j,i+1);
                seen(data(j,i)) = 1;
            else
                data(j,i+1) = dp(data(j,i));
                seen(data(j,i)) = 1;
            end
        end
    end
    disp('x');
end