function F_array = system_equations(x)
    f1 = exp(x(1)-x(2)) - sin(x(1)+x(2));
    f2 = (x(1)*x(2))^2 - cos(x(1)+x(2));
    F_array = [f1; f2];
end