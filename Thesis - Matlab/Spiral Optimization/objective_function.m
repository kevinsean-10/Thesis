function res = objective_function(x)
    f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
    f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
    F_array = [f1,f2];
    res = 0;
    for i = 1:length(F_array)
        res = res + abs(F_array(i));
    end
    res = -1 / (1 + res);
end

