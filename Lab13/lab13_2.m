clear; clc; close all;

% --- ПАРАМЕТРЫ СЕТКИ И ГРАНИЦЫ ---
delta_h = 0.1; % Величина шага интегрирования, заданная в условии
boundary_x1 = 1; % Конечное значение x для первого задания
boundary_x2 = 3; % Конечное значение x для второго задания (исследование жесткости)
nodes_t1 = 0:delta_h:boundary_x1; % Создание сетки узлов для первого задания
nodes_t2 = 0:delta_h:boundary_x2; % Создание сетки узлов для второго задания
num_pts1 = length(nodes_t1); % Количество расчетных точек в первом задании
num_pts2 = length(nodes_t2); % Количество расчетных точек во втором задании

% --- ОПИСАНИЕ МАТЕМАТИЧЕСКИХ МОДЕЛЕЙ (СИСТЕМ УРАВНЕНИЙ) ---
% Используем анонимные функции @(t, u1, u2) для описания правых частей систем
ode_f1_low = @(t, u1, u2) u1*exp(-t^2) + t*u2; % Уравнение y1' для устойчивой системы (Задание 1)
ode_f1_high = @(t, u1, u2) u1*exp(t^2) + t*u2;  % Уравнение y1' для жесткой системы (Задание 2)
ode_f2_base = @(t, u1, u2) 3*t - u1 + 2*u2;    % Уравнение y2', общее для обеих систем

% =========================================================
% БОНУС: АПРИОРНАЯ СТРАТЕГИЯ ВЫБОРА ШАГА
% =========================================================
fprintf('--- АПРИОРНЫЙ АНАЛИЗ УСТОЙЧИВОСТИ ---\n');
t_intervals = 0:0.1:3; % Тестовый интервал для оценки Якобиана
max_abs_eig = 0; % Переменная для хранения максимального собственного числа
for m = 1:length(t_intervals) % Цикл по точкам интервала
    % Формируем матрицу Якобиана J = [df1/dy1, df1/dy2; df2/dy1, df2/dy2]
    Matrix_J = [exp(-t_intervals(m)^2), t_intervals(m); -1, 2]; 
    current_eigs = abs(eig(Matrix_J)); % Вычисляем модули собственных чисел матрицы
    max_abs_eig = max(max_abs_eig, max(current_eigs)); % Ищем глобальный максимум спектрального радиуса
end
h_stability_limit = 2 / max_abs_eig; % Расчет границы устойчивости для явных методов
fprintf('  Максимальное |lambda| в системе = %.2f\n', max_abs_eig);
fprintf('  Критический шаг устойчивости: h_max = %.4f\n', h_stability_limit);
fprintf('  Используемый шаг %.1f (анализ завершен)\n\n', delta_h);

% =========================================================
% ПЕРВЫЙ ЭТАП: РЕШЕНИЕ СИСТЕМЫ 1 (ДО X = 1)
% =========================================================

% 1. ПРЯМОЙ (ЯВНЫЙ) МЕТОД ЭЙЛЕРА
out1_euler = zeros(1,num_pts1); out2_euler = zeros(1,num_pts1); % Выделение памяти под результаты
stiffness_vec1 = zeros(1,num_pts1); % Массив для хранения показателей жесткости на каждом шаге
out1_euler(1)=1; out2_euler(1)=1; % Задание начальных условий y1(0)=1, y2(0)=1
for k = 1:num_pts1-1 % Основной цикл расчета
    % Вычисление значений на следующем шаге по формуле: y_new = y_old + h*f(x, y)
    out1_euler(k+1) = out1_euler(k) + delta_h*ode_f1_low(nodes_t1(k), out1_euler(k), out2_euler(k));
    out2_euler(k+1) = out2_euler(k) + delta_h*ode_f2_base(nodes_t1(k), out1_euler(k), out2_euler(k));
    % Оценка жесткости системы на текущем шаге
    Matrix_J = [exp(-nodes_t1(k)^2), nodes_t1(k); -1, 2]; % Матрица Якобиана
    ev_list = abs(eig(Matrix_J)); % Собственные числа
    stiffness_vec1(k) = max(ev_list)/min(ev_list); % Число жесткости S = max|L|/min|L|
end
stiffness_vec1(end) = stiffness_vec1(end-1); % Дублируем последнее значение жесткости для графика

% 2. МОДИФИЦИРОВАННЫЙ МЕТОД ЭЙЛЕРА (ПРЕДИКТОР-КОРРЕКТОР)
out1_mod = zeros(1,num_pts1); out2_mod = zeros(1,num_pts1); % Память под данные
out1_mod(1)=1; out2_mod(1)=1; % Начальные условия
for k = 1:num_pts1-1 % Цикл
    slope1_1 = ode_f1_low(nodes_t1(k), out1_mod(k), out2_mod(k)); % Наклон в начале шага (f1)
    slope1_2 = ode_f2_base(nodes_t1(k), out1_mod(k), out2_mod(k)); % Наклон в начале шага (f2)
    tmp_u1 = out1_mod(k) + delta_h*slope1_1; % Прогноз (предиктор) для y1
    tmp_u2 = out2_mod(k) + delta_h*slope1_2; % Прогноз (предиктор) для y2
    % Коррекция значения с использованием среднего наклона в начале и конце шага
    out1_mod(k+1) = out1_mod(k) + (delta_h/2)*(slope1_1 + ode_f1_low(nodes_t1(k)+delta_h, tmp_u1, tmp_u2));
    out2_mod(k+1) = out2_mod(k) + (delta_h/2)*(slope1_2 + ode_f2_base(nodes_t1(k)+delta_h, tmp_u1, tmp_u2));
end

% 3. КЛАССИЧЕСКИЙ МЕТОД РУНГЕ-КУТТА 4-ГО ПОРЯДКА
out1_rk4 = zeros(1,num_pts1); out2_rk4 = zeros(1,num_pts1); % Память
out1_rk4(1)=1; out2_rk4(1)=1; % Начальные условия
for k = 1:num_pts1-1 % Цикл RK4
    % Вычисление четырех коэффициентов стадий (k1, k2, k3, k4) для системы
    Step_K1 = delta_h * [ode_f1_low(nodes_t1(k), out1_rk4(k), out2_rk4(k)); ode_f2_base(nodes_t1(k), out1_rk4(k), out2_rk4(k))];
    Step_K2 = delta_h * [ode_f1_low(nodes_t1(k)+delta_h/2, out1_rk4(k)+Step_K1(1)/2, out2_rk4(k)+Step_K1(2)/2); ode_f2_base(nodes_t1(k)+delta_h/2, out1_rk4(k)+Step_K1(1)/2, out2_rk4(k)+Step_K1(2)/2)];
    Step_K3 = delta_h * [ode_f1_low(nodes_t1(k)+delta_h/2, out1_rk4(k)+Step_K2(1)/2, out2_rk4(k)+Step_K2(2)/2); ode_f2_base(nodes_t1(k)+delta_h/2, out1_rk4(k)+Step_K2(1)/2, out2_rk4(k)+Step_K2(2)/2)];
    Step_K4 = delta_h * [ode_f1_low(nodes_t1(k)+delta_h, out1_rk4(k)+Step_K3(1), out2_rk4(k)+Step_K3(2)); ode_f2_base(nodes_t1(k)+delta_h, out1_rk4(k)+Step_K3(1), out2_rk4(k)+Step_K3(2))];
    % Финальный расчет значения по взвешенной сумме коэффициентов
    out1_rk4(k+1) = out1_rk4(k) + (Step_K1(1) + 2*Step_K2(1) + 2*Step_K3(1) + Step_K4(1))/6;
    out2_rk4(k+1) = out2_rk4(k) + (Step_K1(2) + 2*Step_K2(2) + 2*Step_K3(2) + Step_K4(2))/6;
end

% 4. КОНТРОЛЬНОЕ РЕШЕНИЕ ВСТРОЕННЫМ ОПЕРАТОРОМ ODE45
[T_ref1, U_ref1] = ode45(@(t,u) [u(1)*exp(-t^2) + t*u(2); 3*t - u(1) + 2*u(2)], nodes_t1, [1;1]);

% ПЕЧАТЬ РЕЗУЛЬТАТОВ В ТЕРМИНАЛ
fprintf('Результаты Задания 1 (x=1):\n');
fprintf('  Явный Эйлер: y1=%.6f, y2=%.6f\n', out1_euler(end), out2_euler(end));
fprintf('  Число жёсткости S(0.5) = %.2f\n\n', stiffness_vec1(6));

% =========================================================
% ВТОРОЙ ЭТАП: РЕШЕНИЕ ЖЕСТКОЙ СИСТЕМЫ (ДО X = 3)
% =========================================================

% 1. ЯВНЫЙ ЭЙЛЕР ДЛЯ ЖЕСТКОЙ ЗАДАЧИ
out1_stiff_eul = zeros(1,num_pts2); out2_stiff_eul = zeros(1,num_pts2);
stiffness_vec2 = zeros(1,num_pts2);
out1_stiff_eul(1)=1; out2_stiff_eul(1)=1;
for k = 1:num_pts2-1
    % Расчет по явному методу (ожидаем потерю устойчивости при h=0.1)
    out1_stiff_eul(k+1) = out1_stiff_eul(k) + delta_h*ode_f1_high(nodes_t2(k), out1_stiff_eul(k), out2_stiff_eul(k));
    out2_stiff_eul(k+1) = out2_stiff_eul(k) + delta_h*ode_f2_base(nodes_t2(k), out1_stiff_eul(k), out2_stiff_eul(k));
    % Оценка жесткости (экспонента теперь положительная, S будет расти)
    Matrix_J2 = [exp(nodes_t2(k)^2), nodes_t2(k); -1, 2];
    ev_list2 = abs(eig(Matrix_J2)); stiffness_vec2(k) = max(ev_list2)/min(ev_list2);
end
stiffness_vec2(end) = stiffness_vec2(end-1);

% 2. НЕЯВНЫЙ МЕТОД ЭЙЛЕРА (СИМВОЛЬНОЕ РЕШЕНИЕ СИСТЕМЫ)
syms Var1_next Var2_next Var1_prev Var2_prev Step_H T_next % Символьные переменные
% Формируем систему уравнений: y_next = y_prev + h * f(x_next, y_next)
equation_sys1 = Var1_next == Var1_prev + Step_H*(Var1_next*exp(T_next^2) + T_next*Var2_next);
equation_sys2 = Var2_next == Var2_prev + Step_H*(3*T_next - Var1_next + 2*Var2_next);
Res_Solver = solve([equation_sys1, equation_sys2], [Var1_next, Var2_next]); % Решаем символьно
% Превращаем символьные формулы в быстрые численные функции
get_u1_next = matlabFunction(Res_Solver.Var1_next, 'Vars', [Var1_prev, Var2_prev, Step_H, T_next]);
get_u2_next = matlabFunction(Res_Solver.Var2_next, 'Vars', [Var1_prev, Var2_prev, Step_H, T_next]);

out1_implicit = zeros(1,num_pts2); out2_implicit = zeros(1,num_pts2); % Память
out1_implicit(1)=1; out2_implicit(1)=1; % Начальные условия
for k = 1:num_pts2-1 % Цикл неявного метода
    % Вычисляем значения по аналитическим формулам, полученным через solve
    out1_implicit(k+1) = get_u1_next(out1_implicit(k), out2_implicit(k), delta_h, nodes_t2(k+1));
    out2_implicit(k+1) = get_u2_next(out1_implicit(k), out2_implicit(k), delta_h, nodes_t2(k+1));
end

% ode45 для эталона во втором задании
[T_ref2, U_ref2] = ode45(@(t,u) [u(1)*exp(t^2) + t*u(2); 3*t - u(1) + 2*u(2)], nodes_t2, [1;1]);

% ПЕЧАТЬ РЕЗУЛЬТАТОВ ЗАДАНИЯ 2
fprintf('Результаты Задания 2 (x=3):\n');
fprintf('  Явный Эйлер (мог разойтись): y1=%.2e\n', out1_stiff_eul(end));
fprintf('  Неявный Эйлер (устойчив):    y1=%.2e\n', out1_implicit(end));
fprintf('  Макс. жесткость системы S = %.2f\n\n', max(stiffness_vec2));

% =========================================================
% ПОСТРОЕНИЕ ГРАФИКОВ
% =========================================================

% ФИГУРА 1: СРАВНЕНИЕ МЕТОДОВ ДЛЯ Y1 (ЗАДАНИЕ 1)
figure('Name', 'Task 1: y1 analysis', 'Color', 'w');
plot(nodes_t1, out1_euler, 'b-o', nodes_t1, out1_mod, 'g-s', nodes_t1, out1_rk4, 'r-^', T_ref1, U_ref1(:,1), 'k--');
xlabel('x'); ylabel('y1(x)'); title('Задание 1: методы для y1');
legend('Эйлер', 'Модиф. Эйлер', 'RK4', 'ode45'); grid on;

% ФИГУРА 2: СРАВНЕНИЕ МЕТОДОВ ДЛЯ Y2 (ЗАДАНИЕ 1)
figure('Name', 'Task 1: y2 analysis', 'Color', 'w');
plot(nodes_t1, out2_euler, 'b-o', nodes_t1, out2_mod, 'g-s', nodes_t1, out2_rk4, 'r-^', T_ref1, U_ref1(:,2), 'k--');
xlabel('x'); ylabel('y2(x)'); title('Задание 1: методы для y2');
legend('Эйлер', 'Модиф. Эйлер', 'RK4', 'ode45'); grid on;

% ФИГУРА 3: ЧИСЛО ЖЕСТКОСТИ ДЛЯ ЗАДАНИЯ 1
figure('Name', 'Task 1: Stiffness', 'Color', 'w');
plot(nodes_t1, stiffness_vec1, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('S(x)'); title('Число жёсткости S(x) (Часть 1)'); grid on;

% ФИГУРА 4: ЯВНЫЙ VS НЕЯВНЫЙ (Y1, ЖЕСТКАЯ СИСТЕМА)
figure('Name', 'Task 2: Euler comparison y1', 'Color', 'w');
plot(nodes_t2, out1_stiff_eul, 'b-o', nodes_t2, out1_implicit, 'm-s', T_ref2, U_ref2(:,1), 'k--');
xlabel('x'); ylabel('y1(x)'); title('Задание 2: Явный vs Неявный (y1)');
legend('Явный Эйлер', 'Неявный Эйлер', 'ode45'); grid on;
xlim([0,3]);
ylim([-2e5,2e5]);

% ФИГУРА 5: ЯВНЫЙ VS НЕЯВНЫЙ (Y2, ЖЕСТКАЯ СИСТЕМА)
figure('Name', 'Task 2: Euler comparison y2', 'Color', 'w');
plot(nodes_t2, out2_stiff_eul, 'b-o', nodes_t2, out2_implicit, 'm-s', T_ref2, U_ref2(:,2), 'k--');
xlabel('x'); ylabel('y2(x)'); title('Задание 2: Явный vs Неявный (y2)');
legend('Явный Эйлер', 'Неявный Эйлер', 'ode45'); grid on;
xlim([0,3]);
ylim([-2e5,2e5]);

% ФИГУРА 6: ЧИСЛО ЖЕСТКОСТИ ДЛЯ ЗАДАНИЯ 2
figure('Name', 'Task 2: Stiffness', 'Color', 'w');
plot(nodes_t2, stiffness_vec2, 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('S(x)'); title('Число жёсткости S(x) (Часть 2)'); grid on;