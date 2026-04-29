clear; clc; close all;

%строим график
figure(1);
h1 = ezplot('y^5 - exp(-x) - x - 3', [-1, 2, -1, 2]);
set(h1, 'Color', 'b', 'LineWidth', 2);
hold on;
h2 = ezplot('y + y^2 - exp(x)', [-1, 2, -1, 2]);
set(h2, 'Color', 'r', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('линии уровня f1=0 (синий), f2=0 (красный)');
grid on; axis equal;
legend('f1=0', 'f2=0');

%начальное приближение
x0 = 1.1; y0 = 1.3;
plot(x0, y0, 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
text(x0+0.08, y0+0.08, 'начальная точка');

%fsolve
fprintf('fsolve:\n');

%система в виде f(x)=0
F_sys = @(X) [X(2)^5 - exp(-X(1)) - X(1) - 3;
              X(2) + X(2)^2 - exp(X(1))];

%настройки решателя
options = optimoptions('fsolve', 'Display', 'iter', ...
                       'TolFun', 1e-12, ...
                       'TolX', 1e-12, ...
                       'MaxIter', 100);
[X_fsolve, fval, exitflag, output] = fsolve(F_sys, [x0; y0], options);

fprintf('\nрезультат: x=%.10f, y=%.10f, итераций=%d\n', ...
        X_fsolve(1), X_fsolve(2), output.iterations);
fprintf('невязка: f1=%.2e, f2=%.2e\n\n', fval(1), fval(2));

%символьное решение
fprintf('символьное решение:\n');

syms x y
eq1 = y^5 == exp(-x) + x + 3;
eq2 = y == exp(x) - y^2;
sol_sym = vpasolve([eq1, eq2], [x, y], [x0, y0]);

x_sym = double(sol_sym.x);
y_sym = double(sol_sym.y);

fprintf('результат: x=%.10f, y=%.10f\n', x_sym, y_sym);
fprintf('проверка: f1=%.2e, f2=%.2e\n\n', ...
    y_sym^5 - exp(-x_sym) - x_sym - 3, ...
    y_sym + y_sym^2 - exp(x_sym));

%метод простых итераций
fprintf('метод простых итераций:\n');

tau = 1.0;
%приведение к виду x = phi1(x,y), y = phi2(x,y)
phi1 = @(x,y) log(y + y^2);
phi2 = @(x,y) (exp(-x) + x + 3)^(1/5);

x_iter = x0; y_iter = y0;
max_iter = 200; tol = 1e-10;
iter_count = 0; convergence_ok = true;

fprintf('итерация |        x        |        y        |   ||j||  |   невязка\n');
fprintf('---------|-----------------|-----------------|---------|------------\n');

for k = 1:max_iter
     x_std = phi1(x_iter, y_iter);
     y_std = phi2(x_iter, y_iter);
    
    % релаксация
    x_new = (1 - tau) * x_iter + tau * x_std;
    y_new = (1 - tau) * y_iter + tau * y_std;
    
    %вычисление якобиана для проверки сходимости
    eps_j = 1e-6;
    J11 = (phi1(x_iter+eps_j, y_iter) - phi1(x_iter, y_iter)) / eps_j;
    J12 = (phi1(x_iter, y_iter+eps_j) - phi1(x_iter, y_iter)) / eps_j;
    J21 = (phi2(x_iter+eps_j, y_iter) - phi2(x_iter, y_iter)) / eps_j;
    J22 = (phi2(x_iter, y_iter+eps_j) - phi2(x_iter, y_iter)) / eps_j;
    J_norm = max(abs(eig([J11 J12; J21 J22])));
    
    %невязка
    F1_val = y_iter^5 - exp(-x_iter) - x_iter - 3;
    F2_val = y_iter + y_iter^2 - exp(x_iter);
    residual = norm([F1_val; F2_val]);
    
    fprintf('%9d | %.13f | %.13f | %7.4f | %.2e\n', ...
            k, x_new, y_new, J_norm, residual);
    
    %проверка условия сходимости
    if J_norm >= 1
        fprintf('||j|| = %.4f >= 1 (условие сходимости не выполняется)\n', J_norm);
        convergence_ok = false;
    end
    
    %проверка остановки
    if max(abs([x_new - x_iter, y_new - y_iter])) < tol && residual < tol
        iter_count = k;
        break;
    end
    
    x_iter = x_new;
    y_iter = y_new;
    iter_count = k;
end

fprintf('\nрезультат: x=%.10f, y=%.10f\n', x_iter, y_iter);
fprintf('итераций: %d\n', iter_count);
if convergence_ok
    fprintf('условие сходимости выполнено на всех итерациях\n\n');
else
    fprintf('условие сходимости нарушено\n\n');
end

%метод зейделя
fprintf('метод зейделя:\n');

phi1 = @(x,y) log(y + y^2);
phi2 = @(x,y) (exp(-x) + x + 3)^(1/5);

x_zeidel = x0; y_zeidel = y0;
max_iter = 200; tol = 1e-10;
iter_zeidel = 0; conv_zeidel_ok = true;

fprintf('итерация |        x        |        y        |   ||j||  |   невязка\n');
fprintf('---------|-----------------|-----------------|---------|------------\n');

for k = 1:max_iter
    %отличие зейделя: при вычислении y используем уже обновлённый x
    x_new = phi1(x_zeidel, y_zeidel);
    y_new = phi2(x_new, y_zeidel);
    
    %вычисление якобиана для проверки сходимости
    eps_j = 1e-6;
    J11 = (phi1(x_zeidel+eps_j, y_zeidel) - phi1(x_zeidel, y_zeidel)) / eps_j;
    J12 = (phi1(x_zeidel, y_zeidel+eps_j) - phi1(x_zeidel, y_zeidel)) / eps_j;
    J21 = (phi2(x_zeidel+eps_j, y_zeidel) - phi2(x_zeidel, y_zeidel)) / eps_j;
    J22 = (phi2(x_zeidel, y_zeidel+eps_j) - phi2(x_zeidel, y_zeidel)) / eps_j;
    J_norm = max(abs(eig([J11 J12; J21 J22])));
    
    %невязка
    F1_val = y_zeidel^5 - exp(-x_zeidel) - x_zeidel - 3;
    F2_val = y_zeidel + y_zeidel^2 - exp(x_zeidel);
    residual = norm([F1_val; F2_val]);
    
    fprintf('%9d | %.13f | %.13f | %7.4f | %.2e\n', ...
            k, x_new, y_new, J_norm, residual);
    
    %проверка условия сходимости
    if J_norm >= 1
        fprintf('||j|| = %.4f >= 1 (условие сходимости не выполняется)\n', J_norm);
        conv_zeidel_ok = false;
    end
    
    %проверка остановки
    if max(abs([x_new - x_zeidel, y_new - y_zeidel])) < tol && residual < tol
        iter_zeidel = k;
        break;
    end
    
    x_zeidel = x_new;
    y_zeidel = y_new;
    iter_zeidel = k;
end

fprintf('\nрезультат метода зейделя:\n');
fprintf('x = %.10f, y = %.10f\n', x_zeidel, y_zeidel);
fprintf('итераций: %d\n', iter_zeidel);
if conv_zeidel_ok
    fprintf('условие сходимости выполнено на всех итерациях\n\n');
else
    fprintf('условие сходимости нарушено\n\n');
end

%метод ньютона
fprintf('метод ньютона:\n');

F1_N = @(x,y) y^5 - exp(-x) - x - 3;
F2_N = @(x,y) y + y^2 - exp(x);

%аналитический якобиан
J11 = @(x,y) exp(-x) - 1;
J12 = @(x,y) 5*y^4;
J21 = @(x,y) -exp(x);
J22 = @(x,y) 1 + 2*y;

x_newt = x0; y_newt = y0;
max_iter = 100;
iter_newt = 0;
lambda = 1;

fprintf('итерация |        x        |        y        |   невязка  |    det(j)\n');
fprintf('---------|-----------------|-----------------|------------|------------\n');

for k = 1:max_iter
    F = [F1_N(x_newt, y_newt); F2_N(x_newt, y_newt)];
    residual = norm(F);
    
    J = [J11(x_newt, y_newt), J12(x_newt, y_newt);
         J21(x_newt, y_newt), J22(x_newt, y_newt)];
    det_J = det(J);
    
    fprintf('%9d | %.13f | %.13f | %.2e | %.4e\n', ...
            k, x_newt, y_newt, residual, det_J);
    
    %проверка остановки
    if residual < 1e-12
        iter_newt = k;
        break;
    end
    
    %проверка условия сходимости: якобиан не должен быть вырожден
    if abs(det_J) < 1e-8
        fprintf('якобиан близок к вырожденному, используется метод с релаксацией\n');
        dx = -J \ F;
        lambda = 0.5;
        x_newt = x_newt + lambda * dx(1);
        y_newt = y_newt + lambda * dx(2);
    else
        dx = -J \ F;
        x_newt = x_newt + dx(1);
        y_newt = y_newt + dx(2);
    end
    
    iter_newt = k;
end

fprintf('\nрезультат метода ньютона:\n');
fprintf('x = %.10f, y = %.10f\n', x_newt, y_newt);
fprintf('итераций: %d\n', iter_newt);
fprintf('невязка: %.2e\n\n', norm([F1_N(x_newt, y_newt); F2_N(x_newt, y_newt)]));

%финальный график
x_root = X_fsolve(1);
y_root = X_fsolve(2);

figure(2);
h1 = ezplot('y^5 - exp(-x) - x - 3', [-0.5, 1.5, 0, 1.5]);
set(h1, 'Color', 'b', 'LineWidth', 2);
hold on;
h2 = ezplot('y + y^2 - exp(x)', [-0.5, 1.5, 0, 1.5]);
set(h2, 'Color', 'r', 'LineWidth', 2);

%начальная точка - синяя точка
plot(x0, y0, 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
%корень - красная звезда
plot(x_root, y_root, 'r*', 'MarkerSize', 15, 'LineWidth', 2);

xlabel('x'); ylabel('y');
title('решение системы нелинейных уравнений');
legend('f_1 = 0', 'f_2 = 0', 'начальная точка (синяя)', 'корень (красная *)', 'Location', 'best');
grid on; axis equal;

text(x_root+0.05, y_root+0.05, sprintf('(%.4f, %.4f)', x_root, y_root), ...
     'Color', 'red', 'FontWeight', 'bold');
 
%проверка подстановкой
fprintf('проверка:\n');
fprintf('подстановка корня (%.10f, %.10f) в исходные уравнения:\n', x_root, y_root);
fprintf('y^5 = %.10f\n', y_root^5);
fprintf('exp(-x) + x + 3 = %.10f\n', exp(-x_root) + x_root + 3);
fprintf('разность = %.2e\n\n', y_root^5 - (exp(-x_root) + x_root + 3));

fprintf('y = %.10f\n', y_root);
fprintf('exp(x) - y^2 = %.10f\n', exp(x_root) - y_root^2);
fprintf('разность = %.2e\n', y_root - (exp(x_root) - y_root^2));

%сводная таблица
fprintf('\nсводная таблица:\n');
fprintf('%-25s | %-12s | %-12s | %-10s\n', 'метод', 'x', 'y', 'итерации');
fprintf('--------------------------|--------------|--------------|------------\n');
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'fsolve', X_fsolve(1), X_fsolve(2), output.iterations);
fprintf('%-25s | %-12.8f | %-12.8f | %-10s\n', 'vpasolve', x_sym, y_sym, '-');
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'простые итерации', x_iter, y_iter, iter_count);
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'метод зейделя', x_zeidel, y_zeidel, iter_zeidel);
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'метод ньютона', x_newt, y_newt, iter_newt);