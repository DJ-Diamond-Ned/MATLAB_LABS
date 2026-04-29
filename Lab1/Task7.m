%a)
syms x y
y = x^3 - x
ezplot(y, [-4 4])

%b)
y = sin(1/x^2)
ezplot(y, [-2 2])

%c)
y = tan(x/2)
ezplot(y, [-pi pi])
%axis input in console with his parametrs(axis([-pi pi -10 10])) and put "Enter"

%d)
Y = exp(x)^((-x^2)/2)
W = x^4-x^2
ezplot(Y, [-1.5 1.5])
hold on
ezplot(W, [-1.5 1.5])