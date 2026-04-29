clear; clc; close all;

fun_x = @(x) (x - 2).^7 + 3 - cos(x);
dfun_x = @(x) 7*(x - 2).^6 + sin(x);
d2fun_x = @(x) 42*(x - 2).^5 + cos(x);
eps_val = 1e-6;

% --- Проверка условий для метода половинного деления ---
start_a = 0; end_b = 1; % Задание границ интервала [0, 1]
if fun_x(start_a)*fun_x(end_b) < 0
    fprintf('Метод половинного деления: f(a)*f(b) = %.4f < 0 - сходится\n', fun_x(start_a)*fun_x(end_b));
else
    fprintf('Метод половинного деления: f(a)*f(b) = %.4f > 0 - расходится\n', fun_x(start_a)*fun_x(end_b));
end

% --- Проверка условий для метода простых итераций ---
phi_x = @(x) 2 - (3 - cos(x)).^(1/7); % Выраженная функция x = phi(x) для итераций
dphi_x = @(x) -sin(x) ./ (7 * (3 - cos(x)).^(6/7));
x_init = 0.5;
if abs(dphi_x(x_init)) < 1 % Условие сходимости МПИ: модуль производной меньше 1
    fprintf('Метод простых итераций: |phi''(x0)| = %.4f < 1 - сходится\n', abs(dphi_x(x_init)));
else % Если модуль производной больше или равен 1
    fprintf('Метод простых итераций: |phi''(x0)| = %.4f > 1 - расходится\n', abs(dphi_x(x_init)));
end

% --- Проверка условий для метода хорд ---
if fun_x(x_init)*d2fun_x(x_init) > 0 % Проверка условия f(x0)*f''(x0) > 0
    fprintf('Метод хорд: f(x0)*f''''(x0) = %.4f > 0 - сходимость обеспечена\n', fun_x(x_init)*d2fun_x(x_init));
else
    fprintf('Метод хорд: f(x0)*f''''(x0) = %.4f < 0 - возможна расходимость\n', fun_x(x_init)*d2fun_x(x_init));
end

% --- Проверка условий для метода касательных (Ньютона) ---
if fun_x(x_init)*d2fun_x(x_init) > 0 % Проверка условия выбора начальной точки
    fprintf('Метод касательных: f(x0)*f''''(x0) = %.2f > 0 - сходимость обеспечена\n', fun_x(x_init)*d2fun_x(x_init));
else
    fprintf('Метод касательных: f(x0)*f''''(x0) = %.2f < 0 - возможна расходимость\n', fun_x(x_init)*d2fun_x(x_init));
end

% --- Проверка условий для метода секущих ---
if fun_x(x_init)*d2fun_x(x_init) > 0 % Проверка условия выбора начальной точки
    fprintf('Метод секущих: f(x0)*f''''(x0) = %.2f > 0 - сходимость обеспечена\n', fun_x(x_init)*d2fun_x(x_init));
else
    fprintf('Метод секущих: f(x0)*f''''(x0) = %.2f < 0 - возможна расходимость\n', fun_x(x_init)*d2fun_x(x_init));
end
fprintf('\n');

% --- Решение встроенным методом fzero ---
etalon_root = fzero(fun_x, [0, 1]); % Нахождение точного корня для сравнения

% --- Метод половинного деления ---
a_side = 0; b_side = 1; % Локальные переменные границ отрезка
steps_bisec = 0; % Счетчик итераций
for i = 1:100
    c_mid = (a_side + b_side)/2; % Нахождение середины текущего интервала
    steps_bisec = steps_bisec + 1;
    if abs(fun_x(c_mid)) < 1e-3 || (b_side - a_side)/2 < 1e-3 % Проверка условия остановки для локализации
        break;
    end
    if fun_x(a_side)*fun_x(c_mid) < 0 % Если корень в левой половине
        b_side = c_mid; % Сужаем интервал справа
    else % Если корень в правой половине
        a_side = c_mid; % Сужаем интервал слева
    end
end
root_bisec = c_mid;

% --- Метод простых итераций ---
x_val = 0.5; % Начальное приближение
steps_fixed = 0; % Счетчик итераций МПИ
for i = 1:100
    x_next = phi_x(x_val); % Вычисление следующего значения по формуле x = phi(x)
    steps_fixed = steps_fixed + 1;
    if abs(x_next - x_val) < eps_val % Проверка достижения заданной точности
        break;
    end
    x_val = x_next;
end
root_fixed = x_next; % Запись результата МПИ

% --- Реализация метода хорд ---
xh0 = 0; xh1 = 1;
steps_chord = 0;
for i = 1:100
    xh2 = xh1 - fun_x(xh1)*(xh1 - xh0)/(fun_x(xh1) - fun_x(xh0)); % Вычисление координаты пересечения хорды с осью
    steps_chord = steps_chord + 1; % Увеличение счетчика шагов
    if abs(xh2 - xh1) < eps_val % Проверка условия сходимости по разности
        break 
    end
    xh0 = xh1; xh1 = xh2;
end
root_chord = xh2;

% --- Метод касательных ---
xn = 1.0;
steps_newton = 0;
for i = 1:100
    xn_new = xn - fun_x(xn)/dfun_x(xn); % Формула Ньютона: x = x - f(x)/f'(x)
    steps_newton = steps_newton + 1;
    if abs(xn_new - xn) < eps_val % Проверка сходимости
        break
    end
    xn = xn_new;
end
root_newton = xn_new;

% --- Метод секущих ---
s0 = 0; s1 = 1;
steps_secant = 0;
for i = 1:100
    s2 = s1 - fun_x(s1)*(s1 - s0)/(fun_x(s1) - fun_x(s0)); % Вычисление следующего приближения
    steps_secant = steps_secant + 1;
    if abs(s2 - s1) < eps_val % Проверка критерия остановки
        break
    end
    s0 = s1; s1 = s2;
end
root_secant = s2;

% --- Вывод итоговых данных в консоль ---
fprintf('fzero (эталон): %.4f\n', etalon_root);
fprintf('Половинного деления: %.4f, итераций: %d\n', root_bisec, steps_bisec);
fprintf('Простых итераций: %.4f, итераций: %d\n', root_fixed, steps_fixed);
fprintf('Хорд: %.4f, итераций: %d\n', root_chord, steps_chord);
fprintf('Касательных: %.4f, итераций: %d\n', root_newton, steps_newton);
fprintf('Секущих: %.4f, итераций: %d\n', root_secant, steps_secant);

% --- Построение графика ---
plot_x = linspace(-0.2, 1.2, 500);
figure;
plot(plot_x, fun_x(plot_x), 'b-', 'LineWidth', 1.5); hold on; grid on; 
plot([-0.2, 1.2], [0, 0], 'k--');
plot(root_newton, fun_x(root_newton), 'r*', 'MarkerSize', 10)
plot(x_init, fun_x(x_init), 'g*', 'MarkerSize', 10)
xlabel('x');
ylabel('f(x)');
title('График f(x) = (x-2)^7 + 3 - cos(x)');
legend('f(x)', 'Нулевой уровень', 'Найденный корень', 'Начальная точка');