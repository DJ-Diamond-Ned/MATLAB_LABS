clear; clc;

%Исходные данные
x = -3.59; dx = 0.01;
y = 0.467; dy = 0.001;
z = 563.2; dz = 0.1;
f = x*sin(y)+z^(1/3);

%частные прозводные
df_dx = sin(y);
df_dy = x*cos(y);
df_dz = (1/3)*z^(-2/3);

%Абсолютная погрешность
df = abs(df_dx)*dx + abs(df_dy)*dy + abs(df_dz)*dz;
%Относительная погрешность
delta_f = df/abs(f);

fprintf('f = %.6f\n', f);
fprintf('Абсолютная погрешность df = %.6f\n', df);
fprintf('Относительная погрешность -d-f = %.2e\n', delta_f);