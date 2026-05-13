clear; clc; close all;

%% 1. ИСХОДНЫЕ ДАННЫЕ
% t - время (в тексте x), y - координаты маршрута, p - веса (надежность данных)
t = [5, 6, 7, 8, 10, 11, 12, 13, 14];
y = [347.6, 205.1, 173.5, 4.008, 38.54, 92.25, 273.63, 550.10, 650.19];
p = [0.9, 0.8, 0.7, 0.5, 0.9, 1, 0.7, 0.8, 0.3];

% Точка, в которую нужно произвести залп (время x1 = 9)
t_target = 9;

n = length(t);
fprintf('--- РАСЧЕТ ДЛЯ МАЙОРА ГРОМОВА (Вариант 17) ---\n');

%% 2. ОПРЕДЕЛЕНИЕ СТЕПЕНИ ПОЛИНОМА
% Судя по данным (минимум в районе t=8-9), лучше всего подойдет 2-я или 3-я степень
deg = 2; 
fprintf('Выбрана степень полинома для траектории: %d\n\n', deg);

%% 3. МЕТОД 1: Полином Вандермонда (без весов)
fprintf('Метод Вандермонда:\n');
V = zeros(n, deg+1);
for i = 1:n
    for j = 0:deg
        V(i, j+1) = t(i)^j;
    end
end
a_vand = V \ y';
y_res_vand = 0;
for j = 0:deg
    y_res_vand = y_res_vand + a_vand(j+1) * t_target^j;
end
fprintf('Координата для залпа (t=9): %.4f\n\n', y_res_vand);

%% 4. МЕТОД 2: Polyfit (стандартный Matlab)
fprintf('Метод polyfit:\n');
coeffs_poly = polyfit(t, y, deg);
y_res_poly = polyval(coeffs_poly, t_target);
fprintf('Координата для залпа (t=9): %.4f\n\n', y_res_poly);

%% 5. МЕТОД 3: fminsearch (с учетом весов p)
fprintf('Метод fminsearch (с весами):\n');

% 1. Создаем настройки, чтобы увеличить лимит попыток (как просит ошибка)
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000);

% 2. Минимизируем взвешенную сумму квадратов ошибок
weighted_err = @(c) sum(p .* (polyval(c, t) - y).^2);

% 3. ВАЖНО: Вместо нулей даем в качестве начальной точки результат polyfit
% Это гарантирует, что fminsearch стартует уже рядом с правильным ответом
c_start = polyfit(t, y, deg); 

% 4. Запускаем поиск с новыми настройками
coeffs_weighted = fminsearch(weighted_err, c_start, options);

y_res_weighted = polyval(coeffs_weighted, t_target);
fprintf('Координата для залпа (t=9): %.4f\n\n', y_res_weighted);

%% 6. ОЦЕНКА ТОЧНОСТИ (RMSE)
rmse = @(y_calc) sqrt(mean((y_calc - y).^2));
fprintf('Оценка точности (RMSE):\n');
fprintf('Вандермонд: %.4f\n', rmse(polyval(flip(a_vand'), t)));
fprintf('Polyfit:    %.4f\n', rmse(polyval(coeffs_poly, t)));
fprintf('С весами:   %.4f\n\n', rmse(polyval(coeffs_weighted, t)));

%% 7. ГРАФИК
t_plot = linspace(min(t), max(t), 100);
y_plot_fit = polyval(coeffs_weighted, t_plot);

figure('Color', 'w');
hold on; grid on;
% Данные разведки
errorbar(t, y, 1./p*10, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Данные разведки'); 
% Аппроксимация
plot(t_plot, y_plot_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Траектория (модель)');
% Точка удара
plot(t_target, y_res_weighted, 'bs', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Точка залпа (t=9)');

xlabel('Время (t)');
ylabel('Положение колонны (y)');
title('Расчет точки артобстрела');
legend('Location', 'northeast');

% Подпись результата на графике
text(t_target+0.2, y_res_weighted, sprintf('Y = %.2f', y_res_weighted), 'Color', 'blue', 'FontWeight', 'bold');

fprintf('РЕКОМЕНДАЦИЯ: Огонь по координате y = %.2f\n', y_res_weighted);