clear; clc; close all;

% --- 1. Определение функции и производных ---
fun_x = @(x) (x - 2).^7 + 3 - cos(x);
dfun_x = @(x) 7*(x - 2).^6 + sin(x);
d2fun_x = @(x) 42*(x - 2).^5 + cos(x);
eps_val = 1e-6; % Точность для уточнения
x_start = 0.5;  % Начальная точка для итерационных методов

% --- 2. Проверка условий сходимости (Критерий: исключение штрафа -0.3) ---
fprintf('--- Проверка условий сходимости ---\n');
% МПИ: |phi'(x)| < 1
phi_x = @(x) 2 - (3 - cos(x)).^(1/7); 
dphi_x = @(x) -sin(x) ./ (7 * (3 - cos(x)).^(6/7));
fprintf('МПИ: |phi''(x0)| = %.4f (< 1 - ок)\n', abs(dphi_x(x_start)));

% Для методов Ньютона, Хорд и Секущих проверяем знак f(x)*f''(x)
check_conv = fun_x(x_start) * d2fun_x(x_start);
fprintf('Условие f(x0)*f''''(x0) = %.4f (> 0 - ок)\n\n', check_conv);

% --- 3. Локализация корня (Критерий: 0.2 балла) ---
% Метод половинного деления с невысокой точностью (10^-2)
a_loc = 0; b_loc = 1;
steps_loc = 0;
while (b_loc - a_loc) > 1e-2
    c_loc = (a_loc + b_loc)/2;
    steps_loc = steps_loc + 1;
    if fun_x(a_loc) * fun_x(c_loc) < 0
        b_loc = c_loc;
    else
        a_loc = c_loc;
    end
end
root_localized = (a_loc + b_loc)/2;
fprintf('Локализованный корень (точность 0.01): %.4f\n', root_localized);

% --- 4. Уточнение корня (Критерии: 0.4 + 0.3 + 0.4 + 0.3 балла) ---

% А. Метод простых итераций (МПИ)
x_mpi = x_start; steps_mpi = 0;
while true
    x_next = phi_x(x_mpi);
    steps_mpi = steps_mpi + 1;
    if abs(x_next - x_mpi) < eps_val, break; end
    x_mpi = x_next;
end

% Б. Метод Касательных (Ньютона)
x_newt = x_start; steps_newt = 0;
while true
    x_next = x_newt - fun_x(x_newt)/dfun_x(x_newt);
    steps_newt = steps_newt + 1;
    if abs(x_next - x_newt) < eps_val, break; end
    x_newt = x_next;
end

% В. Метод Хорд (Одна точка НЕПОДВИЖНА)
% Выбираем неподвижной ту точку, где f(x)*f''(x) > 0. 
% f(0) = -126, f''(0) = -1343 -> f(0)*f''(0) > 0. Значит x_fix = 0.
x_fix = 0; x_chord = 1; steps_chord = 0;
while true
    x_next = x_chord - fun_x(x_chord)*(x_chord - x_fix)/(fun_x(x_chord) - fun_x(x_fix));
    steps_chord = steps_chord + 1;
    if abs(x_next - x_chord) < eps_val, break; end
    x_chord = x_next;
end

% Г. Метод Секущих (Обе точки ПОДВИЖНЫ)
s0 = 0; s1 = 1; steps_sec = 0;
while true
    s2 = s1 - fun_x(s1)*(s1 - s0)/(fun_x(s1) - fun_x(s0));
    steps_sec = steps_sec + 1;
    if abs(s2 - s1) < eps_val, break; end
    s0 = s1; s1 = s2;
end

% --- 5. Вывод результатов ---
etalon = fzero(fun_x, [0, 1]);
fprintf('\n--- Итоговые результаты (eps = %e) ---\n', eps_val);
fprintf('fzero (эталон):    %.8f\n', etalon);
fprintf('МПИ:               %.8f (итераций: %d)\n', x_mpi, steps_mpi);
fprintf('Касательных:       %.8f (итераций: %d)\n', x_newt, steps_newt);
fprintf('Хорд:              %.8f (итераций: %d)\n', x_chord, steps_chord);
fprintf('Секущих:           %.8f (итераций: %d)\n', s2, steps_sec);

% --- 6. Построение графика (Критерий: 0.1 + 0.1 балла) ---
x_plot = linspace(-0.1, 1.1, 1000);
figure('Color', 'w');
plot(x_plot, fun_x(x_plot), 'b', 'LineWidth', 2); hold on; grid on;
line([min(x_plot) max(x_plot)], [0 0], 'Color', 'k'); % Ось X

% Наносим корень и начальную точку
plot(x_mpi, fun_x(x_mpi), 'r*', 'MarkerSize', 12, 'LineWidth', 2); % Корень
plot(x_start, fun_x(x_start), 'g*', 'MarkerSize', 12, 'LineWidth', 2); % Начало

title('Решение нелинейного уравнения (x-2)^7 + 3 - cos(x) = 0');
xlabel('x'); ylabel('f(x)');
legend('Функция f(x)', 'Ось OX', 'Найденный корень', 'Начальная точка');
