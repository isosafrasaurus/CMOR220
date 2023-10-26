% double
% char
% randperm
% log
% fileread
% strip
% load

Ms = load('letterprob.mat');
M = Ms.letterprob;
mes = fileread('encodedtext.txt');
mes = strip(mes);

y = 1:27;

y(downlow('sghr`kzb`hm`ozqshatkq`hr`rdmrhshud`sn`lhrszjdr'))
y(downlow(mes));
downlowinv(y(downlow(mes)));

f('sghr`kzb`hm`ozqshatkq`hr`rdmrhshud`sn`lhrszjdr', y, M)

function [val] = downlow(text)
    val = double(text) - 95;
end

function [val] = downlowinv(nums)
    val = char(nums + 95);
end

function [val] = f(input, y, M)
    sum = 0;
    for k=1:length(y)-1
        sum = sum + log(M(y(downlow(input(k))),y(downlow(input(k+1)))));
        disp(y(downlow(input(k))) + " " + y(downlow(input(k+1))));
    end
    val = sum;
end
