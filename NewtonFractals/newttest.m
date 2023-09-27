function newttest
    q = [1, 3, 2];
    prettyPrint(q);
end

function prettyPrint(q)
    n = length(q);
    str = '';
    for i = 1:n
        if real(q(i)) ~= 0
            if n-i == 0
                if real(q(i)) < 0
                    str = [str, '-', num2str(abs(real(q(i))))];
                else
                    str = [str, '+', num2str(real(q(i)))];
                end
            else
                if real(q(i)) == -1
                    str = [str, '-x^', num2str(n-i)];
                elseif real(q(i)) == 1
                    str = [str, '+x^', num2str(n-i)];
                elseif real(q(i)) < 0
                    str = [str, '-', num2str(abs(real(q(i)))), 'x^', num2str(n-i)];
                else
                    str = [str, '+', num2str(real(q(i))), 'x^', num2str(n-i)];
                end
            end
        end
        if imag(q(i)) ~= 0
            if n-i == 0
                if imag(q(i)) < 0
                    str = [str, '-', num2str(abs(imag(q(i)))), 'i'];
                else
                    str = [str, '+', num2str(imag(q(i))), 'i'];
                end
            else
                if imag(q(i)) == -1
                    str = [str, '-ix^', num2str(n-i)];
                elseif imag(q(i)) == 1
                    str = [str, '+ix^', num2str(n-i)];
                elseif imag(q(i)) < 0
                    str = [str, '-', num2str(abs(imag(q(i)))), 'ix^', num2str(n-i)];
                else
                    str = [str, '+', num2str(imag(q(i))), 'ix^', num2str(n-i)];
                end
            end
        end
    end
    % Remove leading plus sign.
    if ~isempty(str) && strcmp(str(1:1), '+')
        str(1:1) = [];
    end
    disp(str);
end


