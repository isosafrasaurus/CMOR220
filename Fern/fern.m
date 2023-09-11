function fern
    simplefern(0.3994);

    % reduce size
    % simplefern(0.2);
    % increase size
    % simplefern(0.5);

    advancedfern();
end

function simplefern(cutoff)
    figure()
    hold on

    z = [1 ; 0];
    for i = 1:4000
        plot(z(1), z(2), '.');

        r = rand();
        if r < cutoff
            z = [0.4 -0.3733 ; 0.06 0.6] * z + [0.3533 ; 0];
        else
            z = [-0.8 -0.1867 ; 0.1371 0.80] * z + [1.1 ; 0.1];
        end
    end
    hold off
end

function advancedfern()
    figure()
    hold on

    z = [1 ; 0];
    for i = 1:4000
        plot(z(1), z(2), '.');

        r = rand();
        if r <= 0.01
            z = [0 0 ; 0 0.16] * z;
        elseif r <= 0.76
            z = [0.85 0.04 ; -0.04 0.85] * z + [0 ; 1.6];
        elseif r <= 0.87
            z = [0.2 -0.26 ; 0.23 0.22] * z + [0 ; 1.6];
        else
            z = [-0.15 0.28 ; 0.26 0.24] * z + [0 ; 0.44];
        end
    end
    hold off
end