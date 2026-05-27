clear; clc; close all;

V = 500;         % объем, м^3
k1 = 0,1;        % теплоотдача дна
k2 = 0,57;       % теплоотдача крышки

% Функция теплопотерь Q(d) = 4V/d + pi*d^2*(k1+k2)/4
Q = @(d) 4*V./d + pi*d.^2.*(k1+k2)/4;

% Аналитическое решение
d_analyt = (8*V/(pi*(k1+k2)))^(1/3);
h_analyt = 4*V/(pi*d_analyt^2);
fprintf('Аналитическое решение: d=%.4f м, h=%.4f м, Q=%.4f\n', d_analyt, h_analyt, Q(d_analyt));

% Интервал поиска
a = 1; b = 20;
while Q(b) < Q(b/2), b = b*2;
end

tol = 1e-4;

% Метод золотого сечения
phi = (1+sqrt(5))/2;
resphi = 2-phi;
A = a; B = b;
C = A + resphi*(B-A);
D = B - resphi*(B-A);
iter_gs = 0;

while abs(B-A) > tol
    if Q(C) < Q(D)
        B = D; D = C;
        C = A + resphi*(B-A);
    else
        A = C; C = D;
        D = B - resphi*(B-A);
    end
    iter_gs = iter_gs + 1;
end
d_gs = (A+B)/2;
h_gs = 4*V/(pi*d_gs^2);
fprintf('Золотое сечение: d=%.6f м, h=%.6f м, Q=%.6f, итераций=%d\n', d_gs, h_gs, Q(d_gs), iter_gs);

% Метод парабол
x0 = a; x2 = b; x1 = (x0+x2)/2;
iter_par = 0;

for i = 1:100
    f0 = Q(x0); f1 = Q(x1); f2 = Q(x2);
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

% Метод Ньютона
dQ = @(d) -4*V./d.^2 + pi*d*(k1+k2)/2;
ddQ = @(d) 8*V./d.^3 + pi*(k1+k2)/2;
d_newt = d_analyt * 0.8;
iter_newt = 0;

for i = 1:100
    d_next = d_newt - dQ(d_newt)/ddQ(d_newt);
    if abs(d_next - d_newt) < tol, break; end
    d_newt = d_next;
    iter_newt = i;
end
h_newt = 4*V/(pi*d_newt^2);
fprintf('Метод Ньютона: d=%.6f м, h=%.6f м, Q=%.6f, итераций=%d\n', d_newt, h_newt, Q(d_newt), iter_newt);

% Встроенные ф-ции MATLAB
opts = optimset('Display', 'off');
d_fminbnd = fminbnd(Q, a, b, opts);
h_fminbnd = 4*V/(pi*d_fminbnd^2);
[d_fminsearch, Q_min] = fminsearch(Q, 10, opts);
h_fminsearch = 4*V/(pi*d_fminsearch^2);
fprintf('fminbnd: d=%.6f м, h=%.6f м, Q=%.6f\n', d_fminbnd, h_fminbnd, Q(d_fminbnd));
fprintf('fminsearch: d=%.6f м, h=%.6f м, Q=%.6f\n', d_fminsearch, h_fminsearch, Q_min);

% Графики
figure('Position', [100 100 1000 800]);

% Функция с точками оптимума
subplot(2,2,1);
d_range = linspace(d_analyt*0.3, d_analyt*2.5, 500);
plot(d_range, Q(d_range), 'b-', 'LineWidth', 2); hold on;
plot(d_gs, Q(d_gs), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(d_par, Q(d_par), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(d_newt, Q(d_newt), 'md', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot(d_fminbnd, Q(d_fminbnd), 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('Диаметр d, м'); ylabel('Теплопотери Q(d)');
title('Функция теплопотерь'); legend('Q(d)', 'Золотое сеч.', 'Парабол', 'Ньютона', 'fminbnd');
grid on;

% Шаги метода золотого сечения
subplot(2,2,2);
A = a; B = b; C = A + resphi*(B-A); D = B - resphi*(B-A);
step_vals = [];
while abs(B-A) > tol
    step_vals = [step_vals, (A+B)/2];
    if Q(C) < Q(D), B = D; D = C; C = A + resphi*(B-A);
    else, A = C; C = D; D = B - resphi*(B-A); end
end
plot(1:length(step_vals), step_vals, 'bo-', 'MarkerSize', 6); hold on;
plot(length(step_vals), d_gs, 'r*', 'MarkerSize', 12);
xlabel('Итерация'); ylabel('Текущий диаметр'); title('Сходимость золотого сечения');
grid on;

% Шаги метода Ньютона
subplot(2,2,3);
d_hist = d_analyt * 0.8;
for i = 1:10, d_hist = [d_hist, d_hist(end) - dQ(d_hist(end))/ddQ(d_hist(end))]; end
plot(1:length(d_hist), d_hist, 'mo-', 'MarkerSize', 6); hold on;
plot(length(d_hist), d_newt, 'r*', 'MarkerSize', 12);
xlabel('Итерация'); ylabel('Текущий диаметр'); title('Сходимость метода Ньютона');
grid on;

% Сравнение методов
subplot(2,2,4);
bar([d_gs, d_par, d_newt, d_fminbnd, d_analyt], 'FaceColor', [0.6 0.6 0.9]);
set(gca, 'XTickLabel', {'Зол.сеч.','Параболы','Ньютона','fminbnd','Аналит.'});
ylabel('Оптимальный диаметр, м'); title('Сравнение методов');
grid on;