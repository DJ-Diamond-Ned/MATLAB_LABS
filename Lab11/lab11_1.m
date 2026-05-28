clear; clc; close all;

% --- 1. ИСХОДНЫЕ ДАННЫЕ (Вариант 12) ---
V = 650;                % Объем бака, м^3
k1 = 0.86;              % Теплоотдача дна
k2 = 1.00;              % Теплоотдача крышки
tol = 1e-6;             % Точность

% Функция теплопотерь Q(d) = 4V/d + pi*d^2*(k1+k2)/4
Q = @(d) 4*V./d + pi*d.^2.*(k1+k2)/4;

% Аналитическое решение (для эталона)
d_analyt = (8*V/(pi*(k1+k2)))^(1/3);
h_analyt = 4*V/(pi*d_analyt^2);
fprintf('Аналитическое решение: d=%.4f м, h=%.4f м, Q=%.4f\n', d_analyt, h_analyt, Q(d_analyt));

% Интервал поиска (берем с запасом, чтобы увидеть минимум)
a = 1; b = 50; 

% --- 2. МЕТОД ЗОЛОТОГО СЕЧЕНИЯ ---
phi = (1+sqrt(5))/2;
resphi = 2-phi;
A_gs = a; B_gs = b;
C = A_gs + resphi*(B_gs-A_gs);
D = B_gs - resphi*(B_gs-A_gs);
iter_gs = 0;
gs_history = []; % Для графика сходимости

while abs(B_gs - A_gs) > tol
    gs_history = [gs_history, (A_gs + B_gs)/2];
    if Q(C) < Q(D)
        B_gs = D; D = C; C = A_gs + resphi*(B_gs-A_gs);
    else
        A_gs = C; C = D; D = B_gs - resphi*(B_gs-A_gs);
    end
    iter_gs = iter_gs + 1;
end
d_gs = (A_gs + B_gs)/2;
h_gs = 4*V/(pi*d_gs^2);
fprintf('Золотое сечение: d=%.6f м, h=%.6f м, Q=%.6f, итераций=%d\n', d_gs, h_gs, Q(d_gs), iter_gs);

% --- 3. МЕТОД ПАРАБОЛ ---
x0 = a; x2 = b; x1 = (x0+x2)/2;
iter_par = 0;
for i = 1:100
    f0 = Q(x0); f1 = Q(x1); f2 = Q(x2);
    % Формула аппроксимации
    A_par = ((f2-f0)*(x1-x0) - (f1-f0)*(x2-x0)) / ((x2^2-x0^2)*(x1-x0) - (x1^2-x0^2)*(x2-x0));
    B_par = (f1-f0 - A_par*(x1^2-x0^2))/(x1-x0);
    x_new = -B_par/(2*A_par);
    
    if abs(x_new - x1) < tol, break; end
    if x_new > x1
        if Q(x_new) < f1, x0 = x1; x1 = x_new; else, x2 = x_new; end
    else
        if Q(x_new) < f1, x2 = x1; x1 = x_new; else, x0 = x_new; end
    end
    iter_par = i;
end
d_par = x1;
h_par = 4*V/(pi*d_par^2);
fprintf('Метод парабол: d=%.6f м, h=%.6f м, Q=%.6f, итераций=%d\n', d_par, h_par, Q(d_par), iter_par);

% --- 4. МЕТОД НЬЮТОНА ---
dQ = @(d) -4*V./d.^2 + pi*d*(k1+k2)/2;   % Первая производная
ddQ = @(d) 8*V./d.^3 + pi*(k1+k2)/2;    % Вторая производная
d_newt = d_analyt * 0.8; % Начальная точка
d_hist_newt = [d_newt];
iter_newt = 0;
for i = 1:100
    d_next = d_newt - dQ(d_newt)/ddQ(d_newt);
    d_hist_newt = [d_hist_newt, d_next];
    if abs(d_next - d_newt) < tol, break; end
    d_newt = d_next;
    iter_newt = i;
end
h_newt = 4*V/(pi*d_newt^2);
fprintf('Метод Ньютона: d=%.6f м, h=%.6f м, Q=%.6f, итераций=%d\n', d_newt, h_newt, Q(d_newt), iter_newt);

% --- 5. ВСТРОЕННЫЕ ФУНКЦИИ MATLAB ---
opts = optimset('Display', 'off');
d_fminbnd = fminbnd(Q, a, b, opts);
h_fminbnd = 4*V/(pi*d_fminbnd^2);
[d_fminsearch, Q_min] = fminsearch(Q, (a+b)/2, opts);
h_fminsearch = 4*V/(pi*d_fminsearch^2);
fprintf('fminbnd: d=%.6f м, h=%.6f м, Q=%.6f\n', d_fminbnd, h_fminbnd, Q(d_fminbnd));
fprintf('fminsearch: d=%.6f м, h=%.6f м, Q=%.6f\n', d_fminsearch, h_fminsearch, Q_min);

% --- 6. ГРАФИКИ (Как у друга) ---
figure('Position', [100 100 1000 800], 'Color', 'w');

% 1. Основная функция теплопотерь
subplot(2,2,1);
d_plot = linspace(a+2, b+5, 500);
plot(d_plot, Q(d_plot), 'b-', 'LineWidth', 2); hold on;
plot(d_gs, Q(d_gs), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(d_par, Q(d_par), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(d_newt, Q(d_newt), 'md', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot(d_fminbnd, Q(d_fminbnd), 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
title('Функция теплопотерь'); xlabel('Диаметр d, м'); ylabel('Теплопотери Q(d)');
legend('Q(d)', 'Золотое сеч.', 'Парабол', 'Ньютона', 'fminbnd'); grid on;

% 2. Сходимость золотого сечения
subplot(2,2,2);
plot(1:length(gs_history), gs_history, 'bo-', 'MarkerSize', 5); hold on;
plot(length(gs_history), d_gs, 'r*', 'MarkerSize', 12);
title('Сходимость золотого сечения'); xlabel('Итерация'); ylabel('Текущий диаметр');
grid on;

% 3. Сходимость Ньютона
subplot(2,2,3);
plot(1:length(d_hist_newt), d_hist_newt, 'mo-', 'MarkerSize', 5); hold on;
plot(length(d_hist_newt), d_newt, 'r*', 'MarkerSize', 12);
title('Сходимость метода Ньютона'); xlabel('Итерация'); ylabel('Текущий диаметр');
grid on;

% 4. Столбчатая диаграмма сравнения
subplot(2,2,4);
bar_vals = [d_gs, d_par, d_newt, d_fminbnd, d_analyt];
bar(bar_vals, 'FaceColor', [0.6 0.6 0.9]);
set(gca, 'XTickLabel', {'Зол.сеч.','Параболы','Ньютона','fminbnd','Аналит.'});
title('Сравнение методов'); ylabel('Оптимальный диаметр, м');
grid on;