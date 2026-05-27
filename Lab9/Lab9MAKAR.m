clear; clc; close all;

% =========================================================
% ЗАДАНИЕ 1: ПРОИЗВОДНЫЕ ФУНКЦИИ, ЗАДАННОЙ ТАБЛИЧНО
% =========================================================
fprintf('=== ЗАДАНИЕ 1: АНАЛИЗ ТАБЛИЧНОЙ ФУНКЦИИ (ДАННЫЕ ЛР №7) ===\n');

% Исходный массив данных из предыдущей лабораторной работы
tabular_x_values = [3.65, 3.84, 4.20, 7.83, 10.48, 12.77, 14.81, 17.21, 20.95];
tabular_y_values = [3.373, 3.471, 3.645, 12.90, 78.736, 90.884, 284.174, 336.146, 807.654];

number_of_points = length(tabular_x_values);
polynomial_degree = number_of_points - 1; % Степень полинома для 9 точек равна 8

% Построение интерполяционного полинома с нормализацией (центрирование и масштабирование)
% normalization_params(1) - среднее значение x, normalization_params(2) - стандартное отклонение
[polynomial_coefficients, ~, normalization_params] = polyfit(tabular_x_values, tabular_y_values, polynomial_degree);

% Получение коэффициентов производных полинома
first_derivative_coefficients = polyder(polynomial_coefficients);
second_derivative_coefficients = polyder(first_derivative_coefficients);

% Определение точек между узлами (1 и 2, а также 8 и 9)
target_points_for_analysis = [
    (tabular_x_values(1) + tabular_x_values(2)) / 2, ...
    (tabular_x_values(8) + tabular_x_values(9)) / 2
];

for idx = 1:length(target_points_for_analysis)
    current_x = target_points_for_analysis(idx);
    
    % Перевод x в нормализованную координату для корректного вычисления polyval
    normalized_x = (current_x - normalization_params(1)) / normalization_params(2);
    
    % Вычисление значений производных с учетом масштабирования независимой переменной
    calculated_first_derivative = (1 / normalization_params(2)) * polyval(first_derivative_coefficients, normalized_x);
    calculated_second_derivative = (1 / normalization_params(2)^2) * polyval(second_derivative_coefficients, normalized_x);
    
    % --- Оценка погрешности первой производной по формуле (13.25) ---
    average_step_h = mean(diff(tabular_x_values));
    
    % Для нахождения 9-й конечной разности (n+1) экстраполируем одну точку
    extra_x_point = tabular_x_values(end) + average_step_h;
    extra_y_point = polyval(polynomial_coefficients, (extra_x_point - normalization_params(1))/normalization_params(2));
    extended_y_values = [tabular_y_values, extra_y_point];
    
    % Вычисление конечной разности 9-го порядка
    finite_difference_9th_order = diff(extended_y_values, 9);
    
    % Произведение разностей (omega) для оценки погрешности
    omega_product_value = prod(current_x - tabular_x_values);
    
    % Вычисление теоретической погрешности (аппроксимация производной n+1 порядка разностью)
    error_estimation_derivative = abs((finite_difference_9th_order(1) / (average_step_h^9 * factorial(9))) / factorial(9) * omega_product_value);
    
    fprintf('Точка анализа x = %.4f:\n', current_x);
    fprintf('  Первая производная y'' = %.6f\n', calculated_first_derivative);
    fprintf('  Вторая производная y'''' = %.6f\n', calculated_second_derivative);
    fprintf('  Оценка погрешности y'': %.4e\n\n', error_estimation_derivative);
end


% =========================================================
% ЗАДАНИЕ 2: ПРОИЗВОДНЫЕ НЕПРЕРЫВНОЙ ФУНКЦИИ
% =========================================================
fprintf('=== ЗАДАНИЕ 2: АНАЛИТИЧЕСКАЯ ФУНКЦИЯ ===\n');

target_function = @(x) cos(x.^2 - 3) + x.^3;
evaluation_point_x = 1; 

% Аналитическое вычисление точных значений производных
% y' = -sin(x^2 - 3)*2x + 3x^2
% y'' = -cos(x^2 - 3)*4x^2 - sin(x^2 - 3)*2 + 6x
exact_first_derivative_value = -sin(evaluation_point_x^2 - 3)*2*evaluation_point_x + 3*evaluation_point_x^2;
exact_second_derivative_value = -cos(evaluation_point_x^2 - 3)*4*evaluation_point_x^2 - sin(evaluation_point_x^2 - 3)*2 + 6*evaluation_point_x;

fprintf('Точные значения в x=1:\n  y'' = %.8f\n  y'''' = %.8f\n\n', ...
    exact_first_derivative_value, exact_second_derivative_value);


% Исследование зависимости точности от величины шага дифференцирования
step_h_values = logspace(-1, -9, 50); 
error_1st_deriv_simple = zeros(size(step_h_values));
error_1st_deriv_multipoint = zeros(size(step_h_values));
error_2nd_deriv_simple = zeros(size(step_h_values));
error_2nd_deriv_multipoint = zeros(size(step_h_values));

for i = 1:length(step_h_values)
    current_h = step_h_values(i);
    
    % Вычисление первой производной: центральная разность O(h^2) и 4-точечная схема O(h^4)
    first_deriv_approx_simple = (target_function(evaluation_point_x + current_h) - target_function(evaluation_point_x - current_h)) / (2 * current_h);
    first_deriv_approx_multi  = (-target_function(evaluation_point_x + 2*current_h) + 8*target_function(evaluation_point_x + current_h) ...
                                 - 8*target_function(evaluation_point_x - current_h) + target_function(evaluation_point_x - 2*current_h)) / (12 * current_h);
    
    error_1st_deriv_simple(i) = abs(first_deriv_approx_simple - exact_first_derivative_value);
    error_1st_deriv_multipoint(i) = abs(first_deriv_approx_multi - exact_first_derivative_value);
    
    % Вычисление второй производной: простая схема O(h^2) и 5-точечная схема O(h^4)
    second_deriv_approx_simple = (target_function(evaluation_point_x + current_h) - 2*target_function(evaluation_point_x) + target_function(evaluation_point_x - current_h)) / current_h^2;
    second_deriv_approx_multi  = (-target_function(evaluation_point_x + 2*current_h) + 16*target_function(evaluation_point_x + current_h) ...
                                  - 30*target_function(evaluation_point_x) + 16*target_function(evaluation_point_x - current_h) ...
                                  - target_function(evaluation_point_x - 2*current_h)) / (12 * current_h^2);
    
    error_2nd_deriv_simple(i) = abs(second_deriv_approx_simple - exact_second_derivative_value);
    error_2nd_deriv_multipoint(i) = abs(second_deriv_approx_multi - exact_second_derivative_value);
end

% Демонстрация результатов для фиксированного шага h = 0.01
[~, index_001] = min(abs(step_h_values - 0.01));
fprintf('Результаты конечно-разностных методов при h = 0.01:\n');
fprintf('  y'' (простая схема):   ошибка = %.2e\n', error_1st_deriv_simple(index_001));
fprintf('  y'' (многоточечная):  ошибка = %.2e\n', error_1st_deriv_multipoint(index_001));
fprintf('  y'''' (простая схема):  ошибка = %.2e\n', error_2nd_deriv_simple(index_001));
fprintf('  y'''' (многоточечная): ошибка = %.2e\n\n', error_2nd_deriv_multipoint(index_001));


% =========================================================
% ВИЗУАЛИЗАЦИЯ РЕЗУЛЬТАТОВ (ГРАФИКИ)
% =========================================================
figure('Position', [100, 100, 1300, 450], 'Color', 'w');

% График погрешности первой производной (логарифмический масштаб)
subplot(1,3,1);
loglog(step_h_values, error_1st_deriv_simple, 'b-o', 'MarkerSize', 3); hold on;
loglog(step_h_values, error_1st_deriv_multipoint, 'r-s', 'MarkerSize', 3);
grid on; xlabel('Величина шага h'); ylabel('Абсолютная погрешность');
title('Точность 1-й производной');
legend('Центральная разность O(h^2)', 'Многоточечная схема O(h^4)', 'Location', 'southwest');

% График погрешности второй производной (логарифмический масштаб)
subplot(1,3,2);
loglog(step_h_values, error_2nd_deriv_simple, 'b-o', 'MarkerSize', 3); hold on;
loglog(step_h_values, error_2nd_deriv_multipoint, 'r-s', 'MarkerSize', 3);
grid on; xlabel('Величина шага h'); ylabel('Абсолютная погрешность');
title('Точность 2-й производной');
legend('Простая схема O(h^2)', 'Многоточечная схема O(h^4)', 'Location', 'southwest');


% График интерполяции табличных данных
subplot(1,3,3);
visual_x_range = linspace(min(tabular_x_values), max(tabular_x_values), 300);
visual_y_values = polyval(polynomial_coefficients, (visual_x_range - normalization_params(1)) / normalization_params(2));
plot(visual_x_range, visual_y_values, 'k-', 'LineWidth', 1.2); hold on;
plot(tabular_x_values, tabular_y_values, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
plot(target_points_for_analysis, polyval(polynomial_coefficients, (target_points_for_analysis - normalization_params(1))/normalization_params(2)), 'gs', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; title('Интерполяционный полином');
xlabel('x'); ylabel('y');
legend('Модель полинома', 'Узлы данных', 'Точки вычисления производных');