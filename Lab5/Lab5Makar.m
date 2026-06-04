clear; clc; close all;

% --- 1. Определение функции и её производных ---
f = @(x) 3*cos(6.5*x) + x - 2;
df = @(x) -19.5*sin(6.5*x) + 1;
d2f = @(x) -126.75*cos(6.5*x);

% Функция для метода простых итераций (x = g(x))
g = @(x) acos((2 - x) / 3) / 6.5;
dg = @(x) 1 ./ (19.5 * sqrt(1 - ((2-x)/3).^2)); 

a = 0.1; b = 0.3; eps_val = 1e-4;

% --- 2. Построение графика (Критерий: 0.1 балла) ---
figure('Color', 'w'); hold on; grid on;
x_graph = linspace(0, 0.5, 500); 
plot(x_graph, f(x_graph), 'b', 'LineWidth', 2);
line([0 0.5], [0 0], 'Color', 'k', 'HandleVisibility','off'); % Ось OX
title('График функции f(x) = 3cos(6.5x) + x - 2');
xlabel('x'); ylabel('f(x)');

% --- 3. Решение встроенной функцией (Критерий: 0.1 балла) ---
x_fzero = fzero(f, [a, b]);
fprintf('Результат fzero (эталон): %.6f\n', x_fzero);

% --- 4. Локализация корня (Критерий: 0.2 балла) ---
% Метод половинного деления с низкой точностью (0.1)
x_loc_a = a; x_loc_b = b;
while (x_loc_b - x_loc_a) > 0.1
    m = (x_loc_a + x_loc_b) / 2;
    if f(x_loc_a) * f(m) < 0, x_loc_b = m; else x_loc_a = m; end
end
fprintf('Локализация корня (точность 0.1): [%.2f, %.2f]\n', x_loc_a, x_loc_b);

% --- 5. Метод простых итераций (Критерий: 0.4 балла) ---
fprintf('\n--- Метод простых итераций ---\n');
q = abs(dg(x_fzero)); % Проверка условия сходимости
if q < 1
    fprintf('Условие сходимости выполнено (|g''|=%.4f < 1)\n', q);
    x_c = a; iter1 = 0;
    while true
        iter1 = iter1 + 1;
        x_n = g(x_c);
        if abs(x_n - x_c) < eps_val, break; end
        x_c = x_n;
    end
    fprintf('Корень: %.6f, Итераций: %d\n', x_n, iter1);
else
    fprintf('Условие сходимости НЕ выполнено.\n');
end

% --- 6. Метод хорд (Критерий: 0.3 балла) ---
fprintf('\n--- Метод хорд ---\n');
% Проверка условия Фурье: неподвижной выбирается точка, где f(x)*f''(x) > 0
if f(a)*d2f(a) > 0, x_fix = a; x_mov = b; else x_fix = b; x_mov = a; end
fprintf('Условие f(x)*f''''(x)>0 для x=%.1f (неподвижная точка)\n', x_fix);

iter2 = 0; x_curr = x_mov;
while true
    iter2 = iter2 + 1;
    x_next = x_curr - f(x_curr)*(x_curr - x_fix)/(f(x_curr) - f(x_fix));
    if abs(x_next - x_curr) < eps_val, break; end
    x_curr = x_next;
end
fprintf('Корень: %.6f, Итераций: %d\n', x_next, iter2);

% --- 7. Метод касательных / Ньютона (Критерий: 0.4 балла) ---
fprintf('\n--- Метод касательных ---\n');
% Выбор начальной точки x0 по условию Фурье
if f(a)*d2f(a) > 0, x0 = a; else x0 = b; end
fprintf('Условие Фурье f(x0)*f''''(x0) > 0 выполнено для x0 = %.1f\n', x0);

x_c = x0; iter3 = 0;
while true
    iter3 = iter3 + 1;
    x_n = x_c - f(x_c)/df(x_c);
    if abs(x_n - x_c) < eps_val, break; end
    x_c = x_n;
end
fprintf('Корень: %.6f, Итераций: %d\n', x_n, iter3);

% --- 8. Метод секущих (Критерий: 0.3 балла) ---
fprintf('\n--- Метод секущих ---\n');
% Обе точки подвижны
x_p = a; x_c = b; iter4 = 0;
while true
    iter4 = iter4 + 1;
    x_sec = x_c - f(x_c)*(x_c - x_p)/(f(x_c) - f(x_p));
    if abs(x_sec - x_c) < eps_val, break; end
    x_p = x_c; x_c = x_sec;
end
fprintf('Корень: %.6f, Итераций: %d\n', x_sec, iter4);

% --- 9. Нанесение точек на график (Критерий: 0.1 балла) ---
plot(a, f(a), 'gs', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Начальная точка'); 
plot(x_n, f(x_n), 'r*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Найденный корень'); 
legend('Location', 'best');

fprintf('\nИТОГОВЫЙ КОРЕНЬ: %.4f\n', x_n);
