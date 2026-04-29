clear; clc; close all;

% --- 1. Исходные данные из таблицы (Вариант №17) ---
x = [0, 0.237, 0.417, 0.483, 0.656, 0.777, 0.817, 1.00, 1.15];
y = [10.0, 10.052, 10.279, 10.836, 11.062, 11.884, 12.019, 13.638, 14.654];
n = length(x); % Количество узлов

% Степень полинома для n точек равна n-1
degree = n - 1; 

% Точки, в которых Аркадию нужно найти значения
x_check = [0.2, 0.7];

% --- 2. Построение таблицы конечных разностей ---
% Это нужно для оценки "гладкости" функции и погрешности
dy_table = zeros(n, n);
dy_table(:, 1) = y';
for j = 2:n
    for i = 1:(n - j + 1)
        dy_table(i, j) = dy_table(i+1, j-1) - dy_table(i, j-1);
    end
end

fprintf('Степень интерполирующего полинома: %d\n', degree);
fprintf('\nТаблица конечных разностей (фрагмент):\n');
disp(dy_table);

% --- 3. Канонический полином через определитель Вандермонда ---
% Решаем систему V * a = y, где V - матрица Вандермонда
V = zeros(n, n);
for i = 1:n
    for j = 1:n
        V(i, j) = x(i)^(j-1);
    end
end
coeff_vand = V \ y'; % Коэффициенты полинома: a0, a1, ..., an-1

% --- 4. Интерполяция стандартными средствами MATLAB (polyfit) ---
% Используем polyfit для получения коэффициентов и polyval для вычислений
% mu используется для улучшения численной устойчивости (центрирование и масштабирование)
[coeff_poly, S, mu] = polyfit(x, y, degree);

% Вычисляем значения в искомых точках x1 и x2
y_check_poly = polyval(coeff_poly, x_check, S, mu);

% --- 5. Интерполяция сплайнами ---
% Кубический сплайн обеспечивает более плавную кривую без осцилляций
spline_model = spline(x, y);
y_check_spline = ppval(spline_model, x_check);

% --- 6. Оценка теоретической погрешности ---
% Используем упрощенную оценку через n-ю конечную разность
max_diff_n = abs(dy_table(1, n)); 
factorial_n = factorial(n);
error_estimate = zeros(1, length(x_check));

for k = 1:length(x_check)
    omega = 1;
    for i = 1:n
        omega = omega * (x_check(k) - x(i));
    end
    % Формула: |R(x)| <= |f^(n)(xi) / n! * omega(x)|
    error_estimate(k) = (max_diff_n / factorial_n) * abs(omega);
end

% --- 7. Визуализация результатов ---
x_plot = linspace(min(x), max(x), 1000); % Сетка для плавного графика
y_poly_plot = polyval(coeff_poly, x_plot, S, mu);
y_spline_plot = ppval(spline_model, x_plot);

figure('Color', 'w', 'Name', 'Интерполяция для Аркадия');

% График полинома
subplot(2,1,1);
plot(x_plot, y_poly_plot, 'b-', 'LineWidth', 1.5); hold on;
plot(x, y, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Узлы');
plot(x_check, y_check_poly, 'gs', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Искомые точки');
title(['Полиномиальная интерполяция (степень ', num2str(degree), ')']);
grid on; legend('Полином', 'Узлы', 'Результат Аркадия', 'Location', 'best');

% График сплайна
subplot(2,1,2);
plot(x_plot, y_spline_plot, 'm-', 'LineWidth', 1.5); hold on;
plot(x, y, 'ro', 'MarkerFaceColor', 'r');
plot(x_check, y_check_spline, 'gs', 'MarkerSize', 8, 'LineWidth', 2);
title('Сплайн-интерполяция (кубическая)');
grid on; legend('Сплайн', 'Узлы', 'Результат Аркадия', 'Location', 'best');

% --- 8. Вывод итогов в консоль ---
fprintf('\nРЕЗУЛЬТАТЫ ДЛЯ АРКАДИЯ:\n');
for i = 1:length(x_check)
    fprintf('Точка x%d = %.1f:\n', i, x_check(i));
    fprintf('  - Значение (полином): %.4f\n', y_check_poly(i));
    fprintf('  - Значение (сплайн):  %.4f\n', y_check_spline(i));
    fprintf('  - Оценочная погрешность: %.2e\n', error_estimate(i));
end