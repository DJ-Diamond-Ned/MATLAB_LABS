clear; clc; close all;

%  Проверка, является ли функция решением ДУ
syms y(x) c
Dy = diff(y, x);
ode = Dy == (y/x)^2 + 2*y/x;
y1(x) = c*x^2/(1 - c*x);
res1 = simplify(subs(ode, y, y1));
disp('Проверка y = c*x^2/(1-c*x):');
disp(res1);
y2(x) = -x;
res2 = simplify(subs(ode, y, y2));
disp('Проверка y = -x:');
disp(res2);

%  Решение задачи Коши
f = @(x,y) 2 - y./x;
x0 = 1; y0 = 0; x_end = 1.5;
y_exact = @(x) x - 1./x;

% Параметры шагов
h1 = 0.1; 
h2 = 0.025;

% Метод Эйлера (h=0.1)
x_euler1 = x0:h1:x_end;
y_euler1 = zeros(size(x_euler1));
y_euler1(1) = y0;
for i = 1:length(x_euler1)-1
    y_euler1(i+1) = y_euler1(i) + h1 * f(x_euler1(i), y_euler1(i));
end

% Метод Эйлера (h=0.025)
x_euler2 = x0:h2:x_end;
y_euler2 = zeros(size(x_euler2));
y_euler2(1) = y0;
for i = 1:length(x_euler2)-1
    y_euler2(i+1) = y_euler2(i) + h2 * f(x_euler2(i), y_euler2(i));
end

% Метод Рунге-Кутты 4 (h=0.1)
x_rk1 = x0:h1:x_end;
y_rk1 = zeros(size(x_rk1));
y_rk1(1) = y0;
for i = 1:length(x_rk1)-1
    xi = x_rk1(i);
    yi = y_rk1(i);
    k1 = f(xi, yi);
    k2 = f(xi + h1/2, yi + h1/2*k1);
    k3 = f(xi + h1/2, yi + h1/2*k2);
    k4 = f(xi + h1, yi + h1*k3);
    y_rk1(i+1) = yi + h1/6*(k1 + 2*k2 + 2*k3 + k4);
end

% Метод Рунге-Кутты 4 (h=0.025)
x_rk2 = x0:h2:x_end;
y_rk2 = zeros(size(x_rk2));
y_rk2(1) = y0;
for i = 1:length(x_rk2)-1
    xi = x_rk2(i);
    yi = y_rk2(i);
    k1 = f(xi, yi);
    k2 = f(xi + h2/2, yi + h2/2*k1);
    k3 = f(xi + h2/2, yi + h2/2*k2);
    k4 = f(xi + h2, yi + h2*k3);
    y_rk2(i+1) = yi + h2/6*(k1 + 2*k2 + 2*k3 + k4);
end

% Решение стандартными операторами MATLAB
[x_ode45, y_ode45] = ode45(f, [x0, x_end], y0);

% Погрешности в конце интервала
y_exact_end = y_exact(x_end);
err_euler1_end = abs(y_euler1(end) - y_exact_end);
err_euler2_end = abs(y_euler2(end) - y_exact_end);
err_rk1_end = abs(y_rk1(end) - y_exact_end);
err_rk2_end = abs(y_rk2(end) - y_exact_end);

% Оценка погрешности по Рунге
R_euler = abs(y_euler2(end) - y_euler1(end));
R_rk4 = abs(y_rk2(end) - y_rk1(end)) / (2^4 - 1);

fprintf('Эйлер h=0.1: погрешность в конце = %.6e\n', err_euler1_end);
fprintf('Эйлер h=0.025: погрешность в конце = %.6e\n', err_euler2_end);
fprintf('РК4 h=0.1: погрешность в конце = %.6e\n', err_rk1_end);
fprintf('РК4 h=0.025: погрешность в конце = %.6e\n', err_rk2_end);
fprintf('Оценка погрешности по Рунге (Эйлер): %.6e\n', R_euler);
fprintf('Оценка погрешности по Рунге (РК4): %.6e\n', R_rk4);

% График всех решений
figure;
plot(x_euler1, y_euler1, 'r-o', 'LineWidth', 1, 'MarkerSize', 4); hold on;
plot(x_euler2, y_euler2, 'r--s', 'LineWidth', 1, 'MarkerSize', 3);
plot(x_rk1, y_rk1, 'b-^', 'LineWidth', 1, 'MarkerSize', 4);
plot(x_rk2, y_rk2, 'b--d', 'LineWidth', 1, 'MarkerSize', 3);
plot(x_ode45, y_ode45, 'k-', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('Решения задачи Коши');
legend('Эйлер h=0.1','Эйлер h=0.025','РК4 h=0.1','РК4 h=0.025','ode45','Location','best');
grid on;
hold off;

% График абсолютной погрешности
y_exact_euler1 = y_exact(x_euler1);
y_exact_euler2 = y_exact(x_euler2);
y_exact_rk1 = y_exact(x_rk1);
y_exact_rk2 = y_exact(x_rk2);

figure;
semilogy(x_euler1, abs(y_euler1 - y_exact_euler1), 'r-o', 'LineWidth', 1, 'MarkerSize', 4); hold on;
semilogy(x_euler2, abs(y_euler2 - y_exact_euler2), 'r--s', 'LineWidth', 1, 'MarkerSize', 3);
semilogy(x_rk1, abs(y_rk1 - y_exact_rk1), 'b-^', 'LineWidth', 1, 'MarkerSize', 4);
semilogy(x_rk2, abs(y_rk2 - y_exact_rk2), 'b--d', 'LineWidth', 1, 'MarkerSize', 3);
xlabel('x'); ylabel('Абсолютная погрешность');
title('Поведение абсолютной погрешности');
legend('Эйлер h=0.1','Эйлер h=0.025','РК4 h=0.1','РК4 h=0.025','Location','best');
grid on;
hold off;