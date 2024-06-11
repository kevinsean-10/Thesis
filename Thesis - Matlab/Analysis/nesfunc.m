classdef nesfunc
    methods (Static)
        % Problem 1
        function F_array = system_equations1(x)
            f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
            f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
            F_array = [f1; f2];
        end

        % Problem 2
        function F_array = system_equations2(x)
            f1 = 0.5 * sin(x(1) * x(2)) - 0.25 * x(2) / pi - 0.5 * x(1);
            f2 = (1 - 0.25 / pi) * (exp(2 * x(1)) - exp(1)) + exp(1) * x(2) / pi - 2 * exp(1) * x(1);
            F_array = [f1; f2];
        end

        % Problem 7
        function F_array = system_equations3(x)
            f1 = x(1)^2-x(1)-x(2)^2-x(2)+x(3)^2;
            f2 = sin(x(2)-exp(x(1)));
            f3 = x(3)-log(abs(x(2)));
            F_array = [f1; f2; f3];
        end
    end
end