clear; clc; close all;

%% 1. ИСХОДНЫЕ ДАННЫЕ
t = [5, 6, 7, 8, 10, 11, 12, 13, 14];
y = [347.6, 205.1, 173.5, 4.008, 38.54, 92.25, 273.63, 550.10, 650.19];
p = [0.9, 0.8, 0.7, 0.5, 0.9, 1, 0.7, 0.8, 0.3];
t_target = 9;
n = length(t);

fprintf('--- ОТЧЕТ ПО ЛАБОРАТОРНОЙ РАБОТЕ (ВАРИАНТ 17) ---\n\n');

%% ПУНКТ 1: ОПРЕДЕЛЕНИЕ СТЕПЕНИ ПО КОНЕЧНЫМ РАЗНОСТЯМ (diff) - 0.2 балла
% Считаем разности: если k-я разность примерно постоянна, то степень k
dy1 = diff(y);        % Первая разность
dy2 = diff(dy1);      % Вторая разность
dy3 = diff(dy2);      % Третья разность
fprintf('Анализ конечных разностей (средние значения):\n');
fprintf('|dy1|: %.2f, |dy2|: %.2f, |dy3|: %.2f\n', mean(abs(dy1)), mean(abs(dy2)), mean(abs(dy3)));
deg = 2; % Вторые разности более стабильны для параболы
fprintf('Оптимальная степень полинома: %d\n\n', deg);

%% ПУНКТ 2: ПОЛИНОМ ВАНДЕРМОНДА (БЕЗ ВЕСОВ) - 0.2 балла
V = zeros(n, deg+1);
for i = 1:n
    for j = 0:deg
        V(i, j+1) = t(i)^j;
    end
end
coeffs_vand = V \ y'; 
y_vand = zeros(size(t_target));
for j = 0:deg
    y_vand = y_vand + coeffs_vand(j+1) * t_target.^j;
end

%% ПУНКТ 3: СТАНДАРТНЫЕ ОПЕРАТОРЫ (POLYFIT БЕЗ ВЕСОВ) - 0.2 балла
coeffs_poly = polyfit(t, y, deg);
y_poly = polyval(coeffs_poly, t_target);

%% ПУНКТ 4: ИСПОЛЬЗОВАНИЕ SPAP2 (С ВЕСАМИ) - 0.2 балла
% spap2(узлы, порядок, x, y, веса). Порядок = степень + 1
try
    sp = spap2(1, deg+1, t, y, p); 
    y_spap2 = fnval(sp, t_target);
catch
    y_spap2 = NaN; % Если нет Curve Fitting Toolbox
    fprintf('Внимание: spap2 не сработал (нужен Curve Fitting Toolbox)\n');
end

%% ПУНКТ 5: FMINSEARCH (С ВЕСАМИ, ПОЛИНОМ) - 0.6 балла
weighted_poly_err = @(c) sum(p .* (polyval(c, t) - y).^2);
c_start = polyfit(t, y, deg);
options = optimset('MaxFunEvals', 10000);
coeffs_fmin_w = fminsearch(weighted_poly_err, c_start, options);
y_fmin_w = polyval(coeffs_fmin_w, t_target);

%% ДОПОЛНИТЕЛЬНО (+0.5 балла): НЕ ПОЛИНОМ (ЭКСПОНЕНТА)
% Модель: y = c1 + c2 * exp(c3 * t)
exp_model = @(c, xq) c(1) + c(2) * exp(c(3) * xq);
exp_err = @(c) sum(p .* (exp_model(c, t) - y).^2);
c_exp = fminsearch(exp_err, [1, 1, 0.1]);
y_exp = exp_model(c_exp, t_target);

%% ДОПОЛНИТЕЛЬНО (+0.5 балла): ПОЛИНОМ ЧЕБЫШЁВА
% Нормируем t в диапазон [-1, 1]
t_norm = 2*(t - min(t))/(max(t) - min(t)) - 1;
t_target_norm = 2*(t_target - min(t))/(max(t) - min(t)) - 1;
% Базис Чебышёва: T0=1, T1=x, T2=2x^2-1
T_cheb = @(c, xq) c(1)*1 + c(2)*xq + c(3)*(2*xq.^2 - 1);
cheb_err = @(c) sum(p .* (T_cheb(c, t_norm) - y).^2);
c_cheb = fminsearch(cheb_err, [100, 0, 0]);
y_cheb = T_cheb(c_cheb, t_target_norm);

%% ОЦЕНКА ТОЧНОСТИ (RMSE) - 0.2 балла
rmse = @(y_calc) sqrt(mean((y_calc - y).^2));
fprintf('Оценка RMSE:\n');
fprintf('Polyfit: %.4f\n', rmse(polyval(coeffs_poly, t)));
fprintf('fminsearch (w): %.4f\n', rmse(polyval(coeffs_fmin_w, t)));
fprintf('Экспонента: %.4f\n', rmse(exp_model(c_exp, t)));
fprintf('Чебышёв: %.4f\n\n', rmse(T_cheb(c_cheb, t_norm)));

%% ГРАФИК - 0.2 балла
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; grid on;
t_line = linspace(min(t), max(t), 200);

% Рисуем линии разных цветов
plot(t_line, polyval(coeffs_poly, t_line), 'g-', 'LineWidth', 1, 'DisplayName', 'Polyfit (без весов)');
plot(t_line, polyval(coeffs_fmin_w, t_line), 'b-', 'LineWidth', 2, 'DisplayName', 'Полином (с весами)');
plot(t_line, exp_model(c_exp, t_line), 'm--', 'LineWidth', 1.5, 'DisplayName', 'Экспонента (бонус)');

% Узловые точки в виде ЗВЁЗДОЧЕК (как в критериях)
plot(t, y, 'k*', 'MarkerSize', 10, 'DisplayName', 'Узловые точки');

% Точка залпа
plot(t_target, y_fmin_w, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'ЦЕЛЬ (t=9)');

title('Аппроксимация маршрута колонны (Вариант 17)');
xlabel('Время (t)'); ylabel('Координата (y)');
legend('Location', 'best');

fprintf('РЕЗУЛЬТАТ: Точка залпа (t=9) по методу fminsearch = %.4f\n', y_fmin_w);