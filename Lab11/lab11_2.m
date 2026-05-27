clear; clc; close all;

% Функция: f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
f = @(x,y) (x.^2 + y - 11).^2 + (x + y.^2 - 7).^2;

% Градиент
grad = @(x,y) [4*x*(x^2+y-11) + 2*(x+y^2-7); 2*(x^2+y-11) + 4*y*(x+y^2-7)];

% Параметры
tol = 1e-4;
max_iter = 5000;
alpha = 0.01;

% Начальная точка
x = 4; y = 4;
fprintf('Начальная точка: (%.2f, %.2f), f = %.4f\n', x, y, f(x,y));

% Градиентный спуск
x_hist = x; y_hist = y; f_hist = f(x,y);

for iter = 1:max_iter
    g = grad(x, y);
    if norm(g) < tol, break; end
    
    % Поиск шага
    alpha_temp = alpha;
    for bt = 1:20
        x_new = x - alpha_temp * g(1);
        y_new = y - alpha_temp * g(2);
        if f(x_new, y_new) < f(x, y)
            alpha = alpha_temp * 1.05;
            break;
        else
            alpha_temp = alpha_temp * 0.5;
        end
    end
    
    x = x - alpha * g(1);
    y = y - alpha * g(2);
    
    % Ограничение диапазона [-4.5, 4.5]
    x = max(-4.5, min(4.5, x));
    y = max(-4.5, min(4.5, y));
    
    x_hist = [x_hist, x];
    y_hist = [y_hist, y];
    f_hist = [f_hist, f(x,y)];
end

% Количество итераций (iter - последний номер итерации в цикле)
n_iter = length(f_hist) - 1;  % количество выполненных итераций

fprintf('\nРезультат:\n');
fprintf('Оптимум: (%.8f, %.8f)\n', x, y);
fprintf('Минимум функции: f = %.10f\n', f(x,y));
fprintf('Итераций: %d\n', n_iter);
fprintf('Норма градиента: %.2e\n\n', norm(grad(x,y)));

% Проверка встроенными функциями
f2 = @(v) (v(1)^2 + v(2) - 11)^2 + (v(1) + v(2)^2 - 7)^2;
opts = optimset('Display', 'off');
[x_fms, f_fms] = fminsearch(f2, [4,4], opts);
fprintf('fminsearch: (%.8f, %.8f), f = %.10f\n', x_fms(1), x_fms(2), f_fms);

% Графики
figure('Position', [100 100 1200 500]);

% График 1: Линии уровня и траектория
subplot(1,3,1);
[X, Y] = meshgrid(linspace(1, 4.5, 100), linspace(1, 4.5, 100));
contour(X, Y, f(X,Y), 30, 'LineWidth', 0.8); hold on;
plot(x_hist, y_hist, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot(x_hist(1), y_hist(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x, y, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('x'); ylabel('y'); title('Траектория спуска');
grid on; colorbar;

% График 2: 3D поверхность
subplot(1,3,2);
surf(X, Y, f(X,Y), 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on;
plot3(x_hist, y_hist, f_hist, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
plot3(x, y, f(x,y), 'r*', 'MarkerSize', 10);
xlabel('x'); ylabel('y'); zlabel('f'); title('Поверхность функции');
view(45, 30); grid on;

% График 3: Сходимость
subplot(1,3,3);
semilogy(0:n_iter, f_hist, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3);
xlabel('Итерация'); ylabel('f(x,y)'); title('Сходимость метода');
grid on;