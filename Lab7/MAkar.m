clear; clc; close all;

% ========================= ИСХОДНЫЕ ДАННЫЕ =========================
data_x = [3.65, 3.84, 4.20, 7.83, 10.48, 12.77, 14.81, 17.21, 20.95];
data_y = [3.373, 3.471, 3.645, 12.90, 78.736, 90.884, 284.174, 336.146, 807.654];

% Добавляем точки, в которых нужно предсказать, в начало для построения сплайна
data2_x = [2.2, 2.39, data_x];
data2_y = [23.324, 17.105, data_y];

amount_of_points = length(data_x);   % = 9
x_to_predict = [2.2, 2.39];

% ========================= ТАБЛИЦА КОНЕЧНЫХ РАЗНОСТЕЙ (по исходным 9 точкам) =========================
dy  = diff(data_y);
d2y = diff(dy);
d3y = diff(d2y);
d4y = diff(d3y);
d5y = diff(d4y);
d6y = diff(d5y);
d7y = diff(d6y);
d8y = diff(d7y);

fprintf('Таблица конечных разностей (до 8-го порядка):\n');
fprintf('   i      x(i)      y(i)      dy      d2y      d3y      d4y      d5y      d6y      d7y      d8y\n');
fprintf('-----------------------------------------------------------------------------------------------\n');
for i = 1:amount_of_points
    fprintf('%4d %9.2f %9.3f', i, data_x(i), data_y(i));
    if i <= amount_of_points-1, fprintf(' %8.2f', dy(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-2, fprintf(' %8.2f', d2y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-3, fprintf(' %8.2f', d3y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-4, fprintf(' %8.2f', d4y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-5, fprintf(' %8.2f', d5y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-6, fprintf(' %8.2f', d6y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-7, fprintf(' %8.2f', d7y(i)); else fprintf(' %8s', ''); end
    if i <= amount_of_points-8, fprintf(' %8.2f', d8y(i)); else fprintf(' %8s', ''); end
    fprintf('\n');
end

% ========================= КАНОНИЧЕСКИЙ ПОЛИНОМ (ВАНДЕРМОНД) =========================
V = zeros(amount_of_points, amount_of_points);
for i = 1:amount_of_points
    for j = 1:amount_of_points
        V(i, j) = data_x(i)^(j-1);
    end
end
coeff_vand = V \ data_y';   % не используется дальше, но для отчёта

% ========================= ПОЛИНОМ ЧЕРЕЗ POLYFIT (8-я степень) =========================
[coeff_poly, S, mu] = polyfit(data_x, data_y, amount_of_points-1);
y_poly_pred = polyval(coeff_poly, x_to_predict, S, mu);

% ========================= ОЦЕНКА ПОГРЕШНОСТИ ПОЛИНОМА =========================
h = mean(diff(data_x));                     % средний шаг
max_diff8 = abs(d8y(1));                    % последняя конечная разность
fact_n = factorial(amount_of_points);       % 9!
omega_poly = ones(1, length(x_to_predict));
for k = 1:length(x_to_predict)
    for i = 1:amount_of_points
        omega_poly(k) = omega_poly(k) * (x_to_predict(k) - data_x(i));
    end
end
abs_error_poly = abs(max_diff8 / (h^amount_of_points * fact_n) * omega_poly);
rel_error_poly = abs(abs_error_poly ./ y_poly_pred);

% ========================= ЛИНЕЙНЫЙ СПЛАЙН (order = 2) через spap2 =========================
order = 2;                                  % степень = 1
n_spline = length(data2_x);                 % теперь 11 точек (2.2, 2.39 и исходные 9)
segments = n_spline - 1;                    % для линейной интерполяции нужно n-1 сегментов
spline_model = spap2(order, segments, data2_x, data2_y);
y_spline_pred = fnval(spline_model, x_to_predict);

% ========================= ОЦЕНКА ПОГРЕШНОСТИ ЛИНЕЙНОГО СПЛАЙНА =========================
% (приближённо по исходным 9 точкам)
max_diff2 = abs(d2y(1));                    % вторая конечная разность
fact_2 = factorial(2);                      % 2
omega_spline = ones(1, length(x_to_predict));
for k = 1:length(x_to_predict)
    for i = 1:2
        omega_spline(k) = omega_spline(k) * (x_to_predict(k) - data_x(i));
    end
end
abs_error_spline = abs(max_diff2 / (h^2 * fact_2) * omega_spline);
rel_error_spline = abs(abs_error_spline ./ y_spline_pred);

% ========================= ПОСТРОЕНИЕ ГРАФИКОВ =========================
x_dense_poly = linspace(min(data_x)-1, max(data_x)+1, 500);
y_poly_dense = polyval(coeff_poly, x_dense_poly, S, mu);

% Плотная сетка для сплайна: от самого левого узла (2.2) до самого правого (20.95)
x_dense_spline = linspace(min(data2_x), max(data2_x), 500);
y_spline_dense = fnval(spline_model, x_dense_spline);

figure('Position', [100, 100, 1200, 500], 'Color', 'w');

% График полинома
subplot(1,2,1);
plot(x_dense_poly, y_poly_dense, 'b-', 'LineWidth', 1.5); hold on;
plot(data_x, data_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(x_to_predict, y_poly_pred, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlabel('x'); ylabel('y');
title('Полиномиальная интерполяция (степень 8)');
legend('Интерполяционный полином', 'Узловые точки', 'Заданные точки', 'Location', 'best');
grid on;

% График линейного сплайна (теперь доходит до зелёных точек)
subplot(1,2,2);
plot(x_dense_spline, y_spline_dense, 'm-', 'LineWidth', 1.5); hold on;
plot(data2_x, data2_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(x_to_predict, y_spline_pred, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlabel('x'); ylabel('y');
title('Сплайн-интерполяция (линейный, степень 1)');
legend('Линейный сплайн', 'Узловые точки', 'Заданные точки', 'Location', 'best');
grid on;

% ========================= ВЫВОД РЕЗУЛЬТАТОВ В КОНСОЛЬ =========================
fprintf('\n========== РЕЗУЛЬТАТЫ ИНТЕРПОЛЯЦИИ ==========\n');
fprintf('Точка x     Полином (8 ст.)    Линейный сплайн\n');
for k = 1:length(x_to_predict)
    fprintf('%.4f     %.6f          %.6f\n', x_to_predict(k), y_poly_pred(k), y_spline_pred(k));
end

fprintf('\n========== ОЦЕНКА ПОГРЕШНОСТИ (ПОЛИНОМ) ==========\n');
for k = 1:length(x_to_predict)
    fprintf('x = %.4f: абсолютная = %.4e, относительная = %.4f%%\n', ...
        x_to_predict(k), abs_error_poly(k), rel_error_poly(k)*100);
end

fprintf('\n========== ОЦЕНКА ПОГРЕШНОСТИ (ЛИНЕЙНЫЙ СПЛАЙН) ==========\n');
for k = 1:length(x_to_predict)
    fprintf('x = %.4f: абсолютная = %.4e, относительная = %.4f%%\n', ...
        x_to_predict(k), abs_error_spline(k), rel_error_spline(k)*100);
end

% Итоговое сравнение
fprintf('\n========== ИТОГОВЫЙ АНАЛИЗ ==========\n');
for k = 1:length(x_to_predict)
    if abs_error_poly(k) < abs_error_spline(k)
        fprintf('x = %.4f: Полином теоретически точнее (погрешность %.2e vs %.2e)\n', ...
            x_to_predict(k), abs_error_poly(k), abs_error_spline(k));
    else
        fprintf('x = %.4f: Линейный сплайн теоретически точнее (погрешность %.2e vs %.2e)\n', ...
            x_to_predict(k), abs_error_spline(k), abs_error_poly(k));
    end
end