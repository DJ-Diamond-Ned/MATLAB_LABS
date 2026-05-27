clear; clc; close all;

% =========================================================
% РАЗДЕЛ 1: ИСХОДНЫЕ ДАННЫЕ И ГРАФИК
% =========================================================

syms x % Символьная переменная для аналитических расчетов
target_function_symbolic = cos(x^2 - 3) + x^3; % Исходная функция

% Границы интегрирования
interval_left_border = 0;
interval_right_border = 2; 
integration_interval_length = interval_right_border - interval_left_border;

% Требуемая точность
tolerance_trapezoid = 1e-2; % Точность 0.01 для трапеций
tolerance_simpson = 1e-4;   % Точность 0.0001 для Метода Симпсона

target_function_numerical = matlabFunction(target_function_symbolic);
x_plotting_range = linspace(interval_left_border, interval_right_border, 1000);
figure('Color', 'w', 'Name', 'График исследуемой функции');
plot(x_plotting_range, target_function_numerical(x_plotting_range), 'LineWidth', 2);
grid on; xlabel('x'); ylabel('f(x)');
title('График функции: cos(x^2 - 3) + x^3');


% =========================================================
% РАЗДЕЛ 2: АНАЛИТИЧЕСКОЕ РЕШЕНИЕ (ЭТАЛОН)
% =========================================================
% Рассчёт точного значения интеграла через первообразную (функция int)
analytical_integral_expression = int(target_function_symbolic, x, interval_left_border, interval_right_border);
analytical_integral_value = double(analytical_integral_expression);

fprintf('1. Аналитическое значение (эталонное): %.8f\n\n', analytical_integral_value);

% =========================================================
% РАЗДЕЛ 3: МЕТОД ТРАПЕЦИЙ
% =========================================================
% Для выбора оптимального шага h находим максимум модуля второй производной
second_derivative_symbolic = diff(target_function_symbolic, x, 2);
second_derivative_numerical = matlabFunction(second_derivative_symbolic);

% Находим максимум (M2)
dense_search_grid = linspace(interval_left_border, interval_right_border, 2000);
second_derivative_values = second_derivative_numerical(dense_search_grid);
max_absolute_second_derivative_M2 = max(abs(second_derivative_values));

% Находим шаг h по формуле из критериев оценивания
step_size_trapezoid = sqrt((12 * tolerance_trapezoid) / (integration_interval_length * max_absolute_second_derivative_M2));
% Подсчитываем количество отрезков
number_of_segments_trapezoid = ceil(integration_interval_length / step_size_trapezoid);
% Подогнав под целое число разбиений, находим итоговый размера шага для метода трапеции
step_size_trapezoid_final = integration_interval_length / number_of_segments_trapezoid;

% Вычисляем интеграл через встроенную функцию trapz
nodes_x_trapezoid = interval_left_border : step_size_trapezoid_final : interval_right_border;
nodes_y_trapezoid = target_function_numerical(nodes_x_trapezoid);
integral_value_trapezoid = trapz(nodes_x_trapezoid, nodes_y_trapezoid);

% $ % ПРОЦЕДУРА РУНГЕ % $ %
% Для уточнения считаем интеграл с шагом равным половине шага в методе трапеций ( получается 2*n отрезков)
number_of_segments_doubled = 2 * number_of_segments_trapezoid;
nodes_x_half_step = linspace(interval_left_border, interval_right_border, number_of_segments_doubled + 1);
integral_value_half_step = trapz(nodes_x_half_step, target_function_numerical(nodes_x_half_step));

% Уточняем результат по формуле Рунге
integral_value_runge_corrected = integral_value_half_step + (integral_value_half_step - integral_value_trapezoid) / 3;

fprintf('2. Результаты метода Трапеций (точность 10^-2):\n');
fprintf('   Расчетный шаг h = %.5f (n = %d), M2 = %.4f\n', ...
    step_size_trapezoid_final, number_of_segments_trapezoid, max_absolute_second_derivative_M2);
fprintf('   Базовое значение: %.6f\n', integral_value_trapezoid);
fprintf('   Уточненное по Рунге: %.6f\n', integral_value_runge_corrected);
fprintf('   Фактическая погрешность (Рунге): %.2e\n\n', abs(integral_value_runge_corrected - analytical_integral_value));


% =========================================================
% РАЗДЕЛ 4: ИССЛЕДОВАНИЕ ВЛИЯНИЯ ШАГА h
% =========================================================
step_size_study_array = logspace(-4, -0.5, 20);
absolute_error_study_array = zeros(size(step_size_study_array));

% Рассчитываем различный значения ошибки при увеличении шага интегрирования
for i = 1:length(step_size_study_array)
    segments_count_current = ceil(integration_interval_length / step_size_study_array(i));
    nodes_x_current = linspace(interval_left_border, interval_right_border, segments_count_current + 1);
    integral_temp = trapz(nodes_x_current, target_function_numerical(nodes_x_current));
    % Сравниваю каждый результат с аналитическим эталоном
    absolute_error_study_array(i) = abs(integral_temp - analytical_integral_value);
end

figure('Color', 'w', 'Name', 'Зависимость ошибки от шага');
loglog(step_size_study_array, absolute_error_study_array, 'o-r', 'LineWidth', 2);
grid on; xlabel('Величина шага h'); ylabel('Абсолютная ошибка');
title('Влияние шага на погрешность метода трапеций');


% =========================================================
% РАЗДЕЛ 5: МЕТОД СИМПСОНА (МЕТОД ПАРАБОЛ)
% =========================================================
% В Методе Смпосона применяется четвёртая производная и её максимум (M4)
fourth_derivative_symbolic = diff(target_function_symbolic, x, 4);
fourth_derivative_numerical = matlabFunction(fourth_derivative_symbolic);

fourth_derivative_values = fourth_derivative_numerical(dense_search_grid);
max_absolute_fourth_derivative_M4 = max(abs(fourth_derivative_values));

% Рассчитываем шаг h для метода Симпсона по соответствующей формуле
step_size_simpson = ((180 * tolerance_simpson) / (integration_interval_length * max_absolute_fourth_derivative_M4))^(1/4);
number_of_segments_simpson = ceil(integration_interval_length / step_size_simpson);

% В методе Симпсона количество отрезков n должно быть строго четным
if mod(number_of_segments_simpson, 2) ~= 0
    number_of_segments_simpson = number_of_segments_simpson + 1;
end
step_size_simpson_final = integration_interval_length / number_of_segments_simpson;

% Находим значения функции в узлах
nodes_x_simpson = interval_left_border : step_size_simpson_final : interval_right_border;
Y = target_function_numerical(nodes_x_simpson);

% Рассчитываем значения интеграла
integral_value_simpson = (step_size_simpson_final / 3) * (Y(1) + Y(end) + ...
    4 * sum(Y(2:2:end-1)) + 2 * sum(Y(3:2:end-2)));

fprintf('3. Результаты метода Симпсона (точность 10^-4):\n');
fprintf('   Расчетный шаг h = %.5f (n = %d), M4 = %.4f\n', ...
    step_size_simpson_final, number_of_segments_simpson, max_absolute_fourth_derivative_M4);
fprintf('   Значение по Симпсону: %.8f\n', integral_value_simpson);
fprintf('   Фактическая погрешность: %.2e\n\n', abs(integral_value_simpson - analytical_integral_value));

% =========================================================
% РАЗДЕЛ 6: 
% =========================================================

% Используем встроенную функцию MATLAB для перепроверки всех расчетов
integral_value_matlab_builtin = integral(target_function_numerical, interval_left_border, interval_right_border);
fprintf('4. Сравнение со встроенной функцией integral(): %.8f\n', integral_value_matlab_builtin);
fprintf('   Разница со встроенным методом: %.2e\n\n', abs(integral_value_matlab_builtin - analytical_integral_value));

% Рассчитываем значение неопределенного интеграла
fprintf('5. Выражение неопределенного интеграла:\n');
indefinite_integral_expression = int(target_function_symbolic, x);
pretty(indefinite_integral_expression); % Вывожу формулу в красивом виде

% Расчет несобственного интеграла
syms t
improper_integrand = exp(-t^2); % В качестве подынтегральной выбрана функция Гаусса
improper_integral_result = int(improper_integrand, t, 0, inf);

fprintf('\n6. Несобственный интеграл exp(-t^2) от 0 до бесконечности:\n');
fprintf('   Символьный результат: %s\n', char(improper_integral_result));
fprintf('   Числовой результат: %.8f\n', double(improper_integral_result));