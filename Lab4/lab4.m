clc; clear;

% --- ИСХОДНЫЕ ДАННЫЕ ---
A = [5, 7, 8, 6;
     7, 5, 6, 7;
     8, 6, 5, 8;
     6, 7, 8, 5];

b = [548; 560; 611; 499];

% Точность и константы
eps = 1e-4;
max_iter = 1000;

% Базовые характеристики матрицы
fprintf('=== Характеристики матрицы ===\n');
fprintf('Детерминант: %.4f\n', det(A));
fprintf('Ранг: %d\n', rank(A));
fprintf('Число обусловленности: %.4f\n', cond(A));

% Решение через linsolve
x_lin = linsolve(A, b);
fprintf('\n=== 1. Решение linsolve ===\n');
disp(x_lin');

% Разложение матрицы для итерационных методов:
D = diag(diag(A)); % Диагональная матрица
L = tril(A, -1);   % Строго нижняя треугольная
U = triu(A, 1);    % Строго верхняя треугольная
E = eye(size(A));  % Единичная матрица

% 2. Метод простых итераций
% Формула: B = E - tau*A, F = tau*b
fprintf('\n=== 2. Метод простых итераций ===\n');
tau = 0.05; % Параметр релаксации

B_mpi = E - tau * A;
F_mpi = tau * b;

% Проверка условия сходимости ||B|| < 1
norm_B_mpi = norm(B_mpi, inf);
fprintf('Норма матрицы B (МПИ): %.4f\n', norm_B_mpi);
if norm_B_mpi < 1
    fprintf('Условие сходимости ||B|| < 1 выполнено.\n');
else
    fprintf('ВНИМАНИЕ: ||B|| >= 1, метод может расходиться!\n');
end

% Итерационный процесс
u = zeros(4, 1);
iter_mpi = 0; % Счётчик итераций

for k = 1:max_iter
    iter_mpi = iter_mpi + 1;
    u_new = B_mpi * u + F_mpi;
    
    if norm(u_new - u, inf) < eps || isnan(u_new(1))
        break;
    end
    u = u_new;
end

fprintf('Количество итераций: %d\n', iter_mpi);
disp(u');

% 3. Метод Якоби
fprintf('\n=== 3. Метод Якоби ===\n');

B_jac = -D \ (L + U);
F_jac = D \ b;

% Проверка условия сходимости ||B|| < 1
norm_B_jac = norm(B_jac, inf);
fprintf('Норма матрицы B (Якоби): %.4f\n', norm_B_jac);
if norm_B_jac < 1
    fprintf('Условие сходимости ||B|| < 1 выполнено.\n');
else
    fprintf('ВНИМАНИЕ: ||B|| >= 1, метод расходится!\n');
end

% Итерационный процесс
u = zeros(4, 1);
iter_jac = 0; % Счётчик итераций

for k = 1:max_iter
    iter_jac = iter_jac + 1;
    u_new = B_jac * u + F_jac;
    
    if norm(u_new - u, inf) < eps || isnan(u_new(1))
        break;
    end
    u = u_new;
end

fprintf('Количество итераций: %d\n', iter_jac);
disp(u');

% 4. Метод Зейделя
fprintf('\n=== 4. Метод Зейделя ===\n');

B_sei = -(L + D) \ U;
F_sei = (L + D) \ b;

% Проверка условия сходимости ||B|| < 1
norm_B_sei = norm(B_sei, inf);
fprintf('Норма матрицы B (Зейдель): %.4f\n', norm_B_sei);
if norm_B_sei < 1
    fprintf('Условие сходимости ||B|| < 1 выполнено.\n');
else
    fprintf('ВНИМАНИЕ: ||B|| >= 1, метод расходится!\n');
end

% Итерационный процесс
u = zeros(4, 1);
iter_sei = 0; % Счётчик итераций

for k = 1:max_iter
    iter_sei = iter_sei + 1;
    u_new = B_sei * u + F_sei;
    
    if norm(u_new - u, inf) < eps || isnan(u_new(1))
        break;
    end
    u = u_new;
end

fprintf('Количество итераций: %d\n', iter_sei);
disp(u');
