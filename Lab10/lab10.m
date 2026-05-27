clear; clc; close all;

% =========================================================
% РАЗДЕЛ 1: ИСХОДНЫЕ ДАННЫЕ И ГРАФИК
% =========================================================
syms x % Создаю символьную переменную для аналитических расчетов
target_function_symbolic = sin(x - 1) - x^3 * cos(x); % Моя функция

% Определяю границы интегрирования (выбраны [0, 2] для наглядности)
interval_start = 0;
interval_end = 2; 
integration_interval_length = interval_end - interval_start;

% Задаю требуемую по условию точность (эпсилон)
tolerance_trapezoid = 1e-2; % Точность 0.01 для метода трапеций
tolerance_simpson = 1e-4;   % Точность 0.0001 для метода Симпсона

% Генерирую числовую функцию для расчетов в конкретных точках
target_function_numerical = matlabFunction(target_function_symbolic);

% Рисую график, чтобы посмотреть, как "скачет" функция
x_plotting_range = linspace(interval_start, interval_end, 1000);
figure('Color', 'w', 'Name', 'График функции друга');
plot(x_plotting_range, target_function_numerical(x_plotting_range), 'LineWidth', 2);
grid on; xlabel('x'); ylabel('f(x)');
title('Исследуемая функция: sin(x-1) - x^3 * cos(x)');


% =========================================================
% РАЗДЕЛ 2: АНАЛИТИЧЕСКОЕ РЕШЕНИЕ
% =========================================================
% Вычисляем точный интеграл через формулы (функция int)
analytical_integral_expression = int(target_function_symbolic, x, interval_start, interval_end);
analytical_integral_value = double(analytical_integral_expression);

fprintf(' ШАГ 1: АНАЛИТИКА \n');
fprintf('Точное значение интеграла (эталон): %.8f\n\n', analytical_integral_value);


% =========================================================
% РАЗДЕЛ 3: МЕТОД ТРАПЕЦИЙ С УТОЧНЕНИЕМ ПО РУНГЕ
% =========================================================
% Считаем вторую производную, чтобы определить идеальный шаг h
second_derivative_symbolic = diff(target_function_symbolic, x, 2);
second_derivative_numerical = matlabFunction(second_derivative_symbolic);

% Ищем максимум второй производной (M2) на плотной сетке
dense_search_grid = linspace(interval_start, interval_end, 2000);
second_derivative_values = second_derivative_numerical(dense_search_grid);
max_absolute_second_derivative_M2 = max(abs(second_derivative_values));

% Рассчитываем шаг h по формуле, чтобы гарантировать точность
step_size_trapezoid = sqrt((12 * tolerance_trapezoid) / (integration_interval_length * max_absolute_second_derivative_M2));
% Находим количество отрезков n и подгоняю шаг
number_of_segments_trapezoid = ceil(integration_interval_length / step_size_trapezoid);
step_size_trapezoid_final = integration_interval_length / number_of_segments_trapezoid;

% Считаем интеграл стандартным методом трапеций (trapz)
nodes_x_trapezoid = interval_start : step_size_trapezoid_final : interval_end;
nodes_y_trapezoid = target_function_numerical(nodes_x_trapezoid);
integral_value_trapezoid = trapz(nodes_x_trapezoid, nodes_y_trapezoid);

% --- УТОЧНЕНИЕ ПО РУНГЕ ---
% Считаем еще раз с двойным количеством разбиений (шаг h/2)
number_of_segments_doubled = 2 * number_of_segments_trapezoid;
nodes_x_half_step = linspace(interval_start, interval_end, number_of_segments_doubled + 1);
integral_value_half_step = trapz(nodes_x_half_step, target_function_numerical(nodes_x_half_step));

% Применяем формулу Рунге для повышения точности
integral_value_runge_corrected = integral_value_half_step + (integral_value_half_step - integral_value_trapezoid) / 3;

fprintf(' ШАГ 2: МЕТОД ТРАПЕЦИЙ \n');
fprintf('Расчетный шаг h = %.5f (n = %d разбиений)\n', step_size_trapezoid_final, number_of_segments_trapezoid);
fprintf('Базовый результат: %.6f\n', integral_value_trapezoid);
fprintf('Уточненный по Рунге: %.6f\n', integral_value_runge_corrected);
fprintf('Ошибка после Рунге: %.2e\n\n', abs(integral_value_runge_corrected - analytical_integral_value));


% =========================================================
% РАЗДЕЛ 4: АНАЛИЗ ЗАВИСИМОСТИ ОШИБКИ ОТ ШАГА
% =========================================================
% Берем 20 разных вариантов шага h и смотрю, как меняется погрешность
step_size_study_array = logspace(-4, -0.5, 20);
absolute_error_study_array = zeros(size(step_size_study_array));

for i = 1:length(step_size_study_array)
    n_current = ceil(integration_interval_length / step_size_study_array(i));
    nodes_temp = linspace(interval_start, interval_end, n_current + 1);
    integral_temp = trapz(nodes_temp, target_function_numerical(nodes_temp));
    % Записываем разницу с аналитическим эталоном
    absolute_error_study_array(i) = abs(integral_temp - analytical_integral_value);
end

% Рисуем логарифмический график - это лучший способ показать точность метода
figure('Color', 'w', 'Name', 'График погрешности');
loglog(step_size_study_array, absolute_error_study_array, 's-b', 'LineWidth', 2);
grid on; xlabel('Размер шага h'); ylabel('Абсолютная ошибка');
title('Как уменьшение шага снижает ошибку (Трапеции)');


% =========================================================
% РАЗДЕЛ 5: МЕТОД СИМПСОНА (eps = 10^-4)
% =========================================================
% Для этого метода нужна четвертая производная и её максимум (M4)
fourth_derivative_symbolic = diff(target_function_symbolic, x, 4);
fourth_derivative_numerical = matlabFunction(fourth_derivative_symbolic);

fourth_derivative_values = fourth_derivative_numerical(dense_search_grid);
max_absolute_fourth_derivative_M4 = max(abs(fourth_derivative_values));

% Рассчитываем шаг h для Симпсона (тут корень 4-й степени в формуле)
step_size_simpson = ((180 * tolerance_simpson) / (integration_interval_length * max_absolute_fourth_derivative_M4))^(1/4);
number_of_segments_simpson = ceil(integration_interval_length / step_size_simpson);

% ВАЖНО: В методе Симпсона число отрезков n ОБЯЗАТЕЛЬНО должно быть четным
if mod(number_of_segments_simpson, 2) ~= 0
    number_of_segments_simpson = number_of_segments_simpson + 1;
end
step_size_simpson_final = integration_interval_length / number_of_segments_simpson;

% Получаем узлы и значения функции в них
nodes_x_simpson = interval_start : step_size_simpson_final : interval_end;
nodes_y_simpson = target_function_numerical(nodes_x_simpson);

% Реализуем формулу Симпсона вручную (коэффициенты 1-4-2-4-...-1)
integral_value_simpson = (step_size_simpson_final / 3) * (nodes_y_simpson(1) + nodes_y_simpson(end) + ...
    4 * sum(nodes_y_simpson(2:2:end-1)) + 2 * sum(nodes_y_simpson(3:2:end-2)));

fprintf(' ШАГ 3: МЕТОД СИМПСОНА \n');
fprintf('Расчетный шаг h = %.5f (n = %d разбиений)\n', step_size_simpson_final, number_of_segments_simpson);
fprintf('Результат по Симпсону: %.8f\n', integral_value_simpson);
fprintf('Фактическая ошибка:    %.2e\n\n', abs(integral_value_simpson - analytical_integral_value));


% =========================================================
% РАЗДЕЛ 6: ПРОВЕРКА И БОНУСЫ
% =========================================================

% Используем встроенный решатель MATLAB как финальную проверку
integral_value_matlab_builtin = integral(target_function_numerical, interval_start, interval_end);
fprintf(' ШАГ 4: СТАНДАРТНЫЙ FSOLVE (INTEGRAL) \n');
fprintf('Значение через integral(): %.8f\n\n', integral_value_matlab_builtin);

% Считаем неопределенный интеграл
fprintf(' ШАГ 5: НЕОПРЕДЕЛЕННЫЙ ИНТЕГРАЛ \n');
indefinite_integral_expression = int(target_function_symbolic, x);
pretty(indefinite_integral_expression); % Красивый вывод формулы

% Бонусный расчет несобственного интеграла
syms t
bonus_integrand = exp(-t) * sin(t); % Пример сходящегося несобственного интеграла
bonus_result = int(bonus_integrand, t, 0, inf);

fprintf('\n ШАГ 6: НЕСОБСТВЕННЫЙ ИНТЕГРАЛ \n');
fprintf('Интеграл exp(-t)*sin(t) от 0 до бесконечности = %s\n', char(bonus_result));
fprintf('Числовое значение: %.8f\n', double(bonus_result));

% --- ВИЗУАЛИЗАЦИЯ НЕСОБСТВЕННОГО ИНТЕГРАЛА ---
% Создаем сетку значений t от 0 до 10 (дальше функция уже почти ноль)
t_plot_range = linspace(0, 10, 1000);
y_bonus_numerical = exp(-t_plot_range) .* sin(t_plot_range);

figure('Color', 'w', 'Name', 'График несобственного интеграла');
hold on;

% Закрашиваем область под графиком (визуализация интеграла)
area(t_plot_range, y_bonus_numerical, 'FaceColor', [0.8 0.9 1], 'EdgeColor', 'none', 'DisplayName', 'Интеграл (Площадь)');

plot(t_plot_range, y_bonus_numerical, 'b', 'LineWidth', 2, 'DisplayName', 'exp(-t)*sin(t)');

grid on;
xlabel('t'); ylabel('g(t)');
title(['Визуализация несобственного интеграла. Значение: ', num22str(double(bonus_result), '%.4f')]);
legend('Location', 'northeast');

% Добавляем стрелку, показывающую, что график уходит в бесконечность
annotation('arrow', [0.75 0.85], [0.2 0.2], 'Color', 'r', 'LineWidth', 2);
text(8.5, 0.05, 't \rightarrow \infty', 'FontSize', 12, 'Color', 'r');

hold off;