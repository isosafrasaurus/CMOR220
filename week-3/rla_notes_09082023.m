% Driver
function rla_notes_09082023
    p = 10000;
    for i = 1:10
        p = simulate(p);
    end
    disp(p);

    figure
    xaxis = [];
    yaxis = [];
    z = [1;1];
    for i = 1:100
        xaxis = [xaxis z(1)];
        z = ziterate(z);
        yaxis = [yaxis z(2)];
    end
    plot(xaxis, yaxis, 'Color','g','LineStyle',':','MarkerSize',3)

end

function [y] = ziterate(x)
    r = rand()
    if r <= 0.5
        y = (1 / sqrt(2)) * [1 -1 ; 1 1] * x;
    else
        y = [sqrt(3)/2 1/2; -1/2 sqrt(3)/2] * x;
    end
end

function [y] = simulate(x)
    r = rand();
    if (r <= 0.2)
        x = x * 0.96;
    elseif (r <= 0.8)
        x = x * 1.04;
    else
        x = x * 1.1;
    end
    y = x;
end