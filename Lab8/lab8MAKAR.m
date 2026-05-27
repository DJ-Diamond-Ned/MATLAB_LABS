clear; clc; close all;

% $ % ИСХОДНЫЕ ДАННЫЕ % $ %    % $ % ИСХОДНЫЕ ДАННЫЕ % $ %    % $ % ИСХОДНЫЕ ДАННЫЕ % $ %
data_x = [10.00,     10.01,      10.02,      10.03,      10.04,      10.05,      10.06,      10.07,      10.08];
data_y = [1.000,     13.308,     25.722,     38.212,     50.809,     63.504,     79.185,     87.513,     102.034];
data_p = [0.9,       1.0,        0.9,        0.8,        0.7,        0.5,        0.9,        1,          0.7];
x_pred = 10.036;
amount_of_points = length(data_x);

% $ % ОПРЕДЕЛЕНИЕ СТЕПЕНИ ПОЛИНОМА % $ %
% Расчёт конечных разностей для определения лучшей степени полинома
dy1 = diff(data_y); % dy1 = [
                    %   (13.308 - 1); (25.722 - 13.308); (38.212 - 25.722); (50.809 - 38.212); 
                    %   (63.504 - 50.809); (79.185 - 63.504); (87.513 - 79.185); (102.034 - 87.513)
                    % ]
                    
dy2 = diff(dy1);    % dy2 = [
                    %   (12.414 - 12.308); (12.490 - 12.414); (12.597 - 12.490); (12.695 - 12.597); 
                    %   (15.681 - 12.695); (8.328 - 15.681); (14.5210 - 8.328)
                    % ]

dy3 = diff(dy2);    % dy3 = [
                    %   (0.076 - 0.106); (0.107 - 0.076); (0.098 - 0.107); (2.986 - 0.098); 
                    %   (-7.353 - 2.986); (6.193 - (-7.353))
                    % ]
                    
dy4 = diff(dy3);    % dy4 = [
                    %   (0.031 - (-0.030)); (-0.009 - (0.031)); (2.888 - (-0.009)); (-10.339 - 2.888); 
                    %   (13.546 - (-10.339))
                    % ]
                    
dy5 = diff(dy4);    % dy5 = [
                    %   (-0.040 - 0.061); (2.897 - (-0.040)); (-13.227 - 2.897); (23.885 - (-13.227))
                    % ]
                    
dy6 = diff(dy5);    % dy6 = [
                    %   (2.937 - (-0.101)); (-16.124 - 2.937); (37.112 - -16.124)
                    % ] 
                    
dy7 = diff(dy6);    % dy7 = [
                    %   (-19.061 - 3.038); (53.236 - (-19.061))
                    % ]
                    
dy8 = diff(dy7);    % dy8 = [
                    %   (72.297 - (-22.099))
                    % ]

fprintf('Таблица конечных разностей:\n');
fprintf('   i      x(i)      y(i)      dy      d2y      d3y      d4y      d5y      d6y      d7y      d8y\n');
fprintf('--------------------------------------------------------------------------------------------------\n');
for i = 1:amount_of_points
    fprintf('%4d %9.2f %9.3f', i, data_x(i), data_y(i));
    if (i <= amount_of_points - 1)
        fprintf(' %8.2f', dy1(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 2)
        fprintf(' %8.2f', dy2(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 3)
        fprintf(' %8.3f', dy3(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 4)
        fprintf(' %8.3f', dy4(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 5)
        fprintf(' %8.3f', dy5(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 6)
        fprintf(' %8.3f', dy6(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 7)
        fprintf(' %8.3f', dy7(i))
    else
        fprintf(' %8s', ''); end
    
    if (i <= amount_of_points - 8)
        fprintf(' %8.3f', dy8(i))
    else
        fprintf(' %8s', ''); end
    
    fprintf('\n');
end

fprintf('\nСредние значения конечных разностей:\n');
fprintf('1-й порядок (ср.знач): %.4f\n', mean(dy1));
fprintf('2-й порядок (ср.знач): %.4f\n', mean(dy2));
fprintf('3-й порядок (ср.знач): %.4f\n', mean(dy3));
fprintf('4-й порядок (ср.знач): %.4f\n', mean(dy4));
fprintf('5-й порядок (ср.знач): %.4f\n', mean(dy5));
fprintf('6-й порядок (ср.знач): %.4f\n', mean(dy6));
fprintf('7-й порядок (ср.знач): %.4f\n', mean(dy7));
fprintf('8-й порядок (ср.знач): %.4f\n', mean(dy8));

m = 2;
fprintf('Оптимальная степень полинома: %d\n\n', m);

% $ % КАНОНИЧЕСКИЙ АППРОКСИМИРУЮЩИЙ ПОЛИНОМ (ЧЕРЕЗ МАТРИЦУ ВАНДЕРМОНДА) % $ % 
% Строим прямоугольную матрицу размером  [{{количество_исходных_точек}} x {{оптималная_степень_полинома + 1}}] ([9x3])
% Аппроксимирующий полином имеет вид:
% y = a[1] + (a[2] * x) + (a[3] * x^2) + (a[4] * x^3) + (a[5] * x^4) + (a[6] * x^5) + (a[7] * x^6) + (a[8] * x^7) + (a[9] * x^8)
% Матрица Вандермонда (V) = матрица из строк, являющихся этим полиномом, в который подставлены значения исходных точек (x, y).
% Получается 3 строки имеющих вид:
% a[1] + (a[2] * x) + (a[3] * x^2) + (a[4] * x^3) + (a[5] * x^4) + (a[6] * x^5) + (a[7] * x^6) + (a[8] * x^7) + (a[9] * x^8) - y = 0,
% где x и y - числа из таблицы исходных исходных данных, а значения a[i] - это неизвестные коэффицинеты, которые нужно найти, решив СЛАУ

V = zeros(amount_of_points, m+1);
for i = 1:amount_of_points
    for j = 1:m+1   % Перебираем строки матрицы (и точки из исходных данных)    (первая строка => координаты первой точки)
        V(i, j) = data_x(i)^(j-1);  % Заполняем строку[i] матрицу степенями x[i] от 0 до 3 ( это коэффициенты при неизвестных 'a').
    end
end
coeffs_vand = V \ data_y'; % Решаем образовавшуюся СЛАУ встроенным методом исключения Гаусса. Таким образом находим коэффициенты 'a'
Y_pred_vand = (coeffs_vand(1)) + (coeffs_vand(2) * x_pred^(1)) + (coeffs_vand(3) * x_pred^(2));


% $ % ПОЛИНОМ ЧЕРЕЗ СТАНДАРТНЫЕ ОПЕРАТОРЫ (POLYFIT) % $ %
coeffs_poly = polyfit(data_x, data_y, m);
% С помощью polyval вычисляем значения полученного полинома в точках x_pred ( x_pred = 10.036 )
Y_pred_poly = polyval(coeffs_poly, x_pred);

% $ % ПОЛИНОМ ЧЕРЕЗ СПЛАЙН С ВЕСАМИ (SPAP2) % $ %
order = m+1;
sp = spap2(2, order, data_x, data_y, data_p); 
Y_pred_spap = fnval(sp, x_pred);

% $ % ПОЛИНОМ ЧЕРЕЗ FMINSEARCH С ВЕСАМИ % $ %
% Для стабильности, чтобы иксы лежали около нуля и при возведении в квадрат не давали больших значений,
% используем центрирование (сдвиг x и y на 10). Таким образом улучшаем численную обусловленность
x_centered = data_x - 10;
x_pred_centered = x_pred - 10;

fminsearch_func = @(c, xq) c(1)*xq.^2 + c(2)*xq + c(3);    % Вид нашей аппроксимирующей функции
fminsearch_err = @(c) sum(data_p .* (fminsearch_func(c, x_centered) - data_y).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000, 'Display', 'off');
coeffs_fminsearch = fminsearch(fminsearch_err, [0, 0, 0], options);
Y_pred_fminsearch = fminsearch_func(coeffs_fminsearch, x_pred_centered);

% $ % ОЦЕНКА ТОЧНОСТИ % $ %
rmse_calc = @(y_est) sqrt(mean((y_est - data_y).^2)); % Функция для расчета RMSE (корня среднеквадратичной ошибки)
fprintf('Точность аппроксимации (RMSE):\n');
fprintf('Вандермонд:   %.4f\n', rmse_calc(polyval(flip(coeffs_vand), data_x)));
fprintf('Polyfit:      %.4f\n', rmse_calc(polyval(coeffs_poly, data_x)));
fprintf('spap2:        %.4f\n', rmse_calc(fnval(sp, data_x)));
fprintf('Fminsearch:   %.4f\n\n', rmse_calc(fminsearch_func(coeffs_fminsearch, x_centered)));

% $ % ЭКСПОНЕНЦИАЛЬНАЯ АППРОКСИМАЦИЯ % $ %
exp_func = @(c, xq) c(1) * exp(c(2) * xq) + c(3);
err_exponent = @(c) sum(data_p .* (exp_func(c, x_centered) - data_y).^2);
coeffs_exponent = fminsearch(err_exponent, [0, 0, 0], options);
Y_pred_bonus = exp_func(coeffs_exponent, x_pred_centered);

% $ % ПОЛИНОМ ЧЕБЫШЁВА % $ %
x_min = min(data_x); % x_min будет равно -1
x_max = max(data_x); % x+max будет равно +1
x_normalized = 2*(data_x - x_min)/(x_max - x_min) - 1; % Нормируем значения икс (помещаем в диапазон от [-1 до 1])
% Матрица из полиномов Чебышёва
%(первый столбец: единицы ||| второй столбец: иксы ||| третий столбец: 2*(x^2) - 1)
A_chebyshev = [ones(amount_of_points, 1), x_normalized', 2*x_normalized'.^2 - 1]; 
target_chebyshev = @(c) sum(data_p' .* (A_chebyshev*c' - data_y').^2);
coeffs_chebyshev = fminsearch(target_chebyshev, [50, 50, 1], options);
x_pred_normalized = 2*(x_pred - x_min)/(x_max - x_min) - 1;
Y_pred_chebyshev = coeffs_chebyshev(1) + coeffs_chebyshev(2)*x_pred_normalized + coeffs_chebyshev(3)*(2*x_pred_normalized^2 - 1);

% $ % ПОСТРОЕНИЕ ГРАФИКА % $ %
x_grid = linspace(min(data_x), max(data_x), 300);
x_grid_shifted = x_grid - 10; 
figure('Color', 'w', 'Position', [200 200 800 500]);
hold on; grid on;
plot(data_x, data_y, 'k*', 'MarkerSize', 10, 'DisplayName', 'Узловые точки');
plot(x_grid, polyval(coeffs_poly, x_grid), 'b-', 'LineWidth', 1, 'DisplayName', 'Polyfit');
plot(x_grid, fminsearch_func(coeffs_fminsearch, x_grid_shifted), 'r--', 'LineWidth', 2, 'DisplayName', 'fminsearch');
plot(x_grid, exp_func(coeffs_exponent, x_grid_shifted), 'g:', 'LineWidth', 1.5, 'DisplayName', 'Экспонента');
plot(x_pred, Y_pred_fminsearch, 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'pred');
xlabel('X'); ylabel('Y');
title('Аппроксимация функции');
legend('Location', 'northwest');

% $ % ФИНАЛЬНЫЙ ВЫВОД В КОНСОЛЬ % $ %
fprintf('Рассчётные значения игрек в точке предсказания(x = %.3f):\n', x_pred);
fprintf('y (fminsearch с весами): %.4f\n', Y_pred_fminsearch);
fprintf('y (Чебышёв с весами):   %.4f\n', Y_pred_chebyshev);
fprintf('y (Экспонента):         %.4f\n', Y_pred_bonus);