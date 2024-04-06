function F=FObj(y)
    f1=exp(y(1)-y(2))-sin(y(1)+y(2));
    f2=y(1)^2*y(2)^2-cos(y(1)+y(2));
    F=1/(1+abs(f1)+abs(f2));
