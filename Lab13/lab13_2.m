clear; clc; close all;

% --- ОБЩИЕ ПАРАМЕТРЫ ---
step_h = 0.1;
x_end_task1 = 1;
x_end_task2 = 3;
x_grid1 = 0:step_h:x_end_task1;
x_grid2 = 0:step_h:x_end_task2;
n1 = length(x_grid1);
n2 = length(x_grid2);

% Функции системы
f1_stable = @(x, y1, y2) y1*exp(-x^2) + x*y2; % Задание 1
f1_stiff = @(x, y1, y2) y1*exp(x^2) + x*y2;   % Задание 2
f2_common = @(x, y1, y2) 3*x - y1 + 2*y2;     % Общее второе

% =========================================================
% АПРИОРНАЯ СТРАТЕГИЯ ВЫБОРА ШАГА
% =========================================================
fprintf('априорная стратегия выбора шага:\n');
x_sample = 0:0.1:3;
max_lambda = 0;
for i = 1:length(x_sample)
    J = [exp(-x_sample(i)^2), x_sample(i); -1, 2];
    eig_vals = abs(eig(J));
    max_lambda = max(max_lambda, max(eig_vals));
end
h_max_stable = 2 / max_lambda;
fprintf('  максимальное |lambda| = %.2f\n', max_lambda);
fprintf('  максимальный устойчивый шаг: h_max = %.4f\n', h_max_stable);
fprintf('  (по заданию используем h = 0.1, явные методы могут быть неустойчивы)\n\n');

% =========================================================
% РАСЧЕТЫ ЗАДАНИЯ 1 (до x = 1)
% =========================================================

% 1. Явный метод Эйлера
y1_e = zeros(1,n1); y2_e = zeros(1,n1); 
stiffness1 = zeros(1,n1);
y1_e(1)=1; y2_e(1)=1;
for i = 1:n1-1
    y1_e(i+1) = y1_e(i) + step_h*f1_stable(x_grid1(i), y1_e(i), y2_e(i));
    y2_e(i+1) = y2_e(i) + step_h*f2_common(x_grid1(i), y1_e(i), y2_e(i));
    J = [exp(-x_grid1(i)^2), x_grid1(i); -1, 2];
    ev = abs(eig(J)); stiffness1(i) = max(ev)/min(ev);
end
stiffness1(end) = stiffness1(end-1);

% 2. Модифицированный Эйлер
y1_m = zeros(1,n1); y2_m = zeros(1,n1);
y1_m(1)=1; y2_m(1)=1;
for i = 1:n1-1
    k11 = f1_stable(x_grid1(i), y1_m(i), y2_m(i));
    k12 = f2_common(x_grid1(i), y1_m(i), y2_m(i));
    p1 = y1_m(i) + step_h*k11; p2 = y2_m(i) + step_h*k12;
    y1_m(i+1) = y1_m(i) + (step_h/2)*(k11 + f1_stable(x_grid1(i)+step_h, p1, p2));
    y2_m(i+1) = y2_m(i) + (step_h/2)*(k12 + f2_common(x_grid1(i)+step_h, p1, p2));
end

% 3. Рунге-Кутта 4 порядка
y1_rk = zeros(1,n1); y2_rk = zeros(1,n1);
y1_rk(1)=1; y2_rk(1)=1;
for i = 1:n1-1
    K1 = step_h * [f1_stable(x_grid1(i), y1_rk(i), y2_rk(i)); f2_common(x_grid1(i), y1_rk(i), y2_rk(i))];
    K2 = step_h * [f1_stable(x_grid1(i)+step_h/2, y1_rk(i)+K1(1)/2, y2_rk(i)+K1(2)/2); f2_common(x_grid1(i)+step_h/2, y1_rk(i)+K1(1)/2, y2_rk(i)+K1(2)/2)];
    K3 = step_h * [f1_stable(x_grid1(i)+step_h/2, y1_rk(i)+K2(1)/2, y2_rk(i)+K2(2)/2); f2_common(x_grid1(i)+step_h/2, y1_rk(i)+K2(1)/2, y2_rk(i)+K2(2)/2)];
    K4 = step_h * [f1_stable(x_grid1(i)+step_h, y1_rk(i)+K3(1), y2_rk(i)+K3(2)); f2_common(x_grid1(i)+step_h, y1_rk(i)+K3(1), y2_rk(i)+K3(2))];
    y1_rk(i+1) = y1_rk(i) + (K1(1) + 2*K2(1) + 2*K3(1) + K4(1))/6;
    y2_rk(i+1) = y2_rk(i) + (K1(2) + 2*K2(2) + 2*K3(2) + K4(2))/6;
end

% 4. ode45 (Эталон)
[x_ode1, y_ode1_raw] = ode45(@(x,y) [y(1)*exp(-x^2) + x*y(2); 3*x - y(1) + 2*y(2)], x_grid1, [1;1]);
y1_ode1 = y_ode1_raw(:,1)'; y2_ode1 = y_ode1_raw(:,2)';

% --- Вывод Часть 1 ---
fprintf('явный метод эйлера (часть 1):\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_e(1), y2_e(1));
fprintf('  x=0.5: y1=%.6f, y2=%.6f\n', y1_e(6), y2_e(6));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_e(end), y2_e(end));
fprintf('  число жёсткости на шаге x=0.5: %.2f\n\n', stiffness1(6));

fprintf('модифицированный метод эйлера:\n');
fprintf('  x=1:   y1=%.6f, y2=%.6f\n\n', y1_m(end), y2_m(end));

fprintf('метод рунге-кутта 4 порядка:\n');
fprintf('  x=1:   y1=%.6f, y2=%.6f\n\n', y1_rk(end), y2_rk(end));

% =========================================================
% РАСЧЕТЫ ЗАДАНИЯ 2 (до x = 3)
% =========================================================

% 1. Явный метод Эйлера (Задание 2)
y1_e2 = zeros(1,n2); y2_e2 = zeros(1,n2);
stiffness2 = zeros(1,n2);
y1_e2(1)=1; y2_e2(1)=1;
for i = 1:n2-1
    y1_e2(i+1) = y1_e2(i) + step_h*f1_stiff(x_grid2(i), y1_e2(i), y2_e2(i));
    y2_e2(i+1) = y2_e2(i) + step_h*f2_common(x_grid2(i), y1_e2(i), y2_e2(i));
    J = [exp(x_grid2(i)^2), x_grid2(i); -1, 2];
    ev = abs(eig(J)); stiffness2(i) = max(ev)/min(ev);
end
stiffness2(end) = stiffness2(end-1);

% 2. Неявный метод Эйлера (Символьное решение)
syms Y1n Y2n Y1o Y2o H Xn
eq1 = Y1n == Y1o + H*(Y1n*exp(Xn^2) + Xn*Y2n);
eq2 = Y2n == Y2o + H*(3*Xn - Y1n + 2*Y2n);
S = solve([eq1, eq2], [Y1n, Y2n]);
f_y1_next = matlabFunction(S.Y1n, 'Vars', [Y1o, Y2o, H, Xn]);
f_y2_next = matlabFunction(S.Y2n, 'Vars', [Y1o, Y2o, H, Xn]);

y1_imp = zeros(1,n2); y2_imp = zeros(1,n2);
y1_imp(1)=1; y2_imp(1)=1;
for i = 1:n2-1
    y1_imp(i+1) = f_y1_next(y1_imp(i), y2_imp(i), step_h, x_grid2(i+1));
    y2_imp(i+1) = f_y2_next(y1_imp(i), y2_imp(i), step_h, x_grid2(i+1));
end

% 3. ode45 (Эталон для Задания 2)
[x_ode2, y_ode2_raw] = ode45(@(x,y) [y(1)*exp(x^2) + x*y(2); 3*x - y(1) + 2*y(2)], x_grid2, [1;1]);

% --- Вывод Часть 2 ---
fprintf('явный метод эйлера (часть 2):\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_e2(1), y2_e2(1));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_e2(11), y2_e2(11));
fprintf('  x=3:   y1=%.2e, y2=%.2e\n', y1_e2(end), y2_e2(end));
fprintf('  число жёсткости на шаге x=1: %.2f\n\n', stiffness2(11));

fprintf('неявный метод эйлера (символьное решение):\n');
fprintf('  получены аналитические формулы для y1_new и y2_new\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_imp(1), y2_imp(1));
fprintf('  x=3:   y1=%.6f, y2=%.6f\n\n', y1_imp(end), y2_imp(end));

% =========================================================
% ГРАФИКИ (6 ОТДЕЛЬНЫХ ОКОН)
% =========================================================

% Фигура 1: y1(x) Часть 1
figure('Name', 'часть 1 y1(x)', 'Color', 'w');
plot(x_grid1, y1_e, 'b-o', x_grid1, y1_m, 'g-s', x_grid1, y1_rk, 'r-^', x_ode1, y1_ode1, 'k--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y1(x)'); title('часть 1: сравнение методов y1(x)');
legend('явный эйлер', 'модиф. эйлер', 'рунге-кутта 4', 'ode45', 'Location', 'best'); grid on;

% Фигура 2: y2(x) Часть 1
figure('Name', 'часть 1 y2(x)', 'Color', 'w');
plot(x_grid1, y2_e, 'b-o', x_grid1, y2_m, 'g-s', x_grid1, y2_rk, 'r-^', x_ode1, y2_ode1, 'k--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y2(x)'); title('часть 1: сравнение методов y2(x)');
legend('явный эйлер', 'модиф. эйлер', 'рунге-кутта 4', 'ode45', 'Location', 'best'); grid on;

% Фигура 3: Жёсткость Часть 1
figure('Name', 'часть 1 число жёсткости', 'Color', 'w');
plot(x_grid1, stiffness1, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('число жёсткости'); title('число жёсткости на каждом шаге (часть 1)'); grid on;

% Фигура 4: y1(x) Часть 2
figure('Name', 'часть 2 y1(x)', 'Color', 'w');
plot(x_grid2, y1_e2, 'b-o', x_grid2, y1_imp, 'm-s', x_ode2, y_ode2_raw(:,1), 'k--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y1(x)'); title('часть 2: явный vs неявный эйлер y1(x)');
legend('явный эйлер', 'неявный эйлер (символьн.)', 'ode45', 'Location', 'best'); grid on;

% Фигура 5: y2(x) Часть 2
figure('Name', 'часть 2 y2(x)', 'Color', 'w');
plot(x_grid2, y2_e2, 'b-o', x_grid2, y2_imp, 'm-s', x_ode2, y_ode2_raw(:,2), 'k--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y2(x)'); title('часть 2: явный vs неявный эйлер y2(x)');
legend('явный эйлер', 'неявный эйлер (символьн.)', 'ode45', 'Location', 'best'); grid on;

% Фигура 6: Жёсткость Часть 2
figure('Name', 'часть 2 число жёсткости', 'Color', 'w');
plot(x_grid2, stiffness2, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('число жёсткости'); title('число жёсткости на каждом шаге (часть 2)'); grid on;