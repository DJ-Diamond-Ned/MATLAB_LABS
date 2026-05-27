clear; clc; close all;

% 1. ИСХОДНЫЕ ДАННЫЕ
t = [5, 6, 7, 8, 10, 11, 12, 13, 14]; % Время (узловые точки)
y = [347.6, 205.1, 173.5, 4.008, 38.54, 92.25, 273.63, 550.10, 650.19]; % Координаты
p = [0.9, 0.8, 0.7, 0.5, 0.9, 1, 0.7, 0.8, 0.3]; % Весовые коэффициенты (доверие к данным)
t_target = 9; % Точка, в которой нужно вычислить y
n = length(t); % Количество экспериментальных точек

fprintf('ВЫВОД РЕЗУЛЬТАТОВ\n\n');

% ПУНКТ 1: ОПРЕДЕЛЕНИЕ СТЕПЕНИ ПО КОНЕЧНЫМ РАЗНОСТЯМ (оператор diff)
dy1 = diff(y);        % Первая разность (y2-y1, y3-y2...)
dy2 = diff(dy1);      % Вторая разность (изменение скорости изменения)
dy3 = diff(dy2);      % Третья разность
dy4 = diff(dy3);      % Четвертая разность
dy5 = diff(dy4);      % Пятая разность
dy6 = diff(dy5);      % Шестая разность

fprintf('Таблица конечных разностей:\n');
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('  i  |   t   |    y     |   dy1    |   dy2    |   dy3    |   dy4    |   dy5    |   dy6    |\n');
fprintf('-------------------------------------------------------------------------------------------\n');

for i = 1:n
    fprintf(' %2d  | %5.1f | %8.3f |', i, t(i), y(i));
    
    if i <= length(dy1)
        fprintf(' %8.3f |', dy1(i));
    else
        fprintf('    -     |');
    end
    
    if i <= length(dy2)
        fprintf(' %8.3f |', dy2(i));
    else
        fprintf('    -     |');
    end
    
    if i <= length(dy3)
        fprintf(' %8.3f |', dy3(i));
    else
        fprintf('    -     |');
    end
    
    if i <= length(dy4)
        fprintf(' %8.3f |', dy4(i));
    else
        fprintf('    -     |');
    end
    
    if i <= length(dy5)
        fprintf(' %8.3f |', dy5(i));
    else
        fprintf('    -     |');
    end
    
    if i <= length(dy6)
        fprintf(' %8.3f |', dy6(i));
    else
        fprintf('    -     |');
    end
    
    fprintf('\n');
end

fprintf('-------------------------------------------------------------------------------------------\n\n');
fprintf('Анализ конечных разностей (средние значения абс. величин):\n');
fprintf('|dy1|: %.2f, |dy2|: %.2f, |dy3|: %.2f\n, |dy4|: %.2f, |dy5|: %.2f, |dy6|: %.2f\n', mean(abs(dy1)), mean(abs(dy2)), mean(abs(dy3)),mean(abs(dy4)), mean(abs(dy5)), mean(abs(dy6)));
deg = 2; % Выбираем 2-ю степень, так как парабола лучше всего описывает разворот траектории
fprintf('Выбранная степень полинома: %d\n\n', deg);

% ПУНКТ 2: ПОЛИНОМ ВАНДЕРМОНДА (БЕЗ ВЕСОВ)
V = zeros(n, deg+1); % Создаем пустую матрицу нужного размера
for i = 1:n
    for j = 0:deg
        V(i, j+1) = t(i)^j; % Заполняем матрицу степенями t (t^0, t^1, t^2)
    end
end
coeffs_vand = V \ y'; % Решаем систему линейных уравнений методом деления матриц
y_vand = 0; % Инициализируем переменную для результата
for j = 0:deg
    y_vand = y_vand + coeffs_vand(j+1) * t_target^j; % Считаем y по формуле полинома
end

% ПУНКТ 3: СТАНДАРТНЫЕ ОПЕРАТОРЫ MATLAB (POLYFIT БЕЗ ВЕСОВ)
coeffs_poly = polyfit(t, y, deg); % Находим коэффициенты стандартной функцией
y_poly = polyval(coeffs_poly, t_target); % Считаем значение в точке 9

% ПУНКТ 4: ИСПОЛЬЗОВАНИЕ SPAP2 (С ВЕСАМИ)
try
    sp = spap2(1, deg+1, t, y, p); % 1 - один сегмент (полином), deg+1 - порядок
    y_spap2 = fnval(sp, t_target); % Вычисляем значение в целевой точке
catch
    y_spap2 = NaN; % Если Toolbox не установлен, пишем "пусто"
    fprintf('Функция spap2 недоступна (нужен Curve Fitting Toolbox)\n');
end

% ПУНКТ 5: FMINSEARCH (С ВЕСАМИ, ПОЛИНОМ)
% Описываем функцию ошибки: сумма весов p, умноженных на квадраты отклонений
weighted_poly_err = @(c) sum(p .* (polyval(c, t) - y).^2);
% Настройки: увеличиваем лимиты итераций для точности
opts = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);
% Стартуем не с нулей, а с результата polyfit (чтобы fminsearch не "улетел")
c_start = polyfit(t, y, deg);
coeffs_fmin_w = fminsearch(weighted_poly_err, c_start, opts); % Ищем минимум ошибки
y_fmin_w = polyval(coeffs_fmin_w, t_target); % Считаем итог

% БОНУС 1: НЕ ПОЛИНОМИАЛЬНАЯ ФУНКЦИЯ (ЭКСПОНЕНТА)
% Модель: y = c1 + c2 * exp(c3 * t)
exp_model = @(c, xq) c(1) + c(2) * exp(c(3) * xq);
% Функция ошибки для экспоненты с учетом весов p
exp_err = @(c) sum(p .* (exp_model(c, t) - y).^2);
% Запуск поиска коэффициентов для экспоненты
c_exp = fminsearch(exp_err, [1, 1, 0.1], opts);
y_exp = exp_model(c_exp, t_target);

% БОНУС 2: ПОЛИНОМ ЧЕБЫШЁВА
% Полиномы Чебышёва требуют, чтобы X был в диапазоне [-1, 1]
t_norm = 2*(t - min(t))/(max(t) - min(t)) - 1; % Нормируем время
t_target_norm = 2*(t_target - min(t))/(max(t) - min(t)) - 1; % Нормируем точку 9
% Описываем базис: c1*T0 + c2*T1 + c3*T2, где T0=1, T1=x, T2=2x^2-1
T_cheb = @(c, xq) c(1)*1 + c(2)*xq + c(3)*(2*xq.^2 - 1);
% Функция ошибки для Чебышёва
cheb_err = @(c) sum(p .* (T_cheb(c, t_norm) - y).^2);
c_cheb = fminsearch(cheb_err, [100, 0, 0], opts); % Поиск коэффициентов
y_cheb = T_cheb(c_cheb, t_target_norm); % Расчет значения

% ОЦЕНКА ТОЧНОСТИ (RMSE)
% Создаем функцию для расчета среднеквадратичной ошибки
rmse_func = @(y_calc) sqrt(mean((y_calc - y).^2));
fprintf('Оценка точности (RMSE):\n');
fprintf('- Polyfit: %.4f\n', rmse_func(polyval(coeffs_poly, t)));
fprintf('- fminsearch (веса): %.4f\n', rmse_func(polyval(coeffs_fmin_w, t)));
fprintf('- Чебышёв (веса): %.4f\n\n', rmse_func(T_cheb(c_cheb, t_norm)));

% ПОСТРОЕНИЕ ГРАФИКА
figure('Color', 'w', 'Name', 'Результаты аппроксимации'); % Белый фон окна
hold on; grid on; % Рисуем всё на одном поле с сеткой
t_line = linspace(min(t), max(t), 200); % Сетка для плавной отрисовки линий

% Рисуем разные модели разными цветами (пункт 7 критериев)
plot(t_line, polyval(coeffs_poly, t_line), 'g:', 'LineWidth', 1.5, 'DisplayName', 'Polyfit (без весов)');
plot(t_line, polyval(coeffs_fmin_w, t_line), 'b-', 'LineWidth', 2, 'DisplayName', 'Полином (с весами)');
plot(t_line, T_cheb(c_cheb, 2*(t_line-min(t))/(max(t)-min(t))-1), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Чебышёв (бонус)');

% Узловые точки в виде ЗВЁЗДОЧЕК (пункт 7 критериев)
plot(t, y, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Узловые точки (звездочки)');

% Точка залпа (наш результат)
plot(t_target, y_fmin_w, 'ms', 'MarkerSize', 12, 'MarkerFaceColor', 'm', 'DisplayName', 'Точка залпа (t=9)');

% Оформление
xlabel('Время (t)'); ylabel('Координата (y)');
title('Аппроксимация траектории для артобстрела');
legend('Location', 'best');

% Вывод финального ответа
fprintf('РЕКОМЕНДАЦИЯ МАЙОРУ: Цель при t=9 находится в точке y = %.4f\n', y_fmin_w);