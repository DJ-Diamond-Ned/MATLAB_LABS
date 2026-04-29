clear; clc; close all;

%строим график
figure(1);
L1 = ezplot('0.1*y^3 - 5*x - y', [-1, 2, -2, 1]);
set(L1, 'Color', 'b', 'LineWidth', 2);
hold on;
L2 = ezplot('y^3 - exp(x) + 3', [-1, 2, -2, 1]);
set(L2, 'Color', 'r', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('линии уровня f1=0 (синий), f2=0 (красный)');
grid on; axis equal;
legend('f1=0', 'f2=0');

%начальное приближение
p0x = 0.6; p0y = -1.1;
plot(p0x, p0y, 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
text(p0x+0.08, p0y+0.08, 'начальная точка');

%fsolve
fprintf('fsolve:\n');

%система в виде f(x)=0
F_system = @(V) [0.1*V(2)^3 - 5*V(1) - V(2);
                 V(2)^3 - exp(V(1)) + 3];

%настройки решателя
opts = optimoptions('fsolve', 'Display', 'iter', ...
                       'TolFun', 1e-12, ...
                       'TolX', 1e-12, ...
                       'MaxIter', 100);
[Res_fsolve, f_vals, exit_f, out_info] = fsolve(F_system, [p0x; p0y], opts);

fprintf('\nрезультат: x=%.10f, y=%.10f, итераций=%d\n', ...
        Res_fsolve(1), Res_fsolve(2), out_info.iterations);
fprintf('невязка: f1=%.2e, f2=%.2e\n\n', f_vals(1), f_vals(2));

%символьное решение
fprintf('символьное решение:\n');

syms sx sy
eqn1 = 0.1*sy^3 - 5*sx - sy == 0;
eqn2 = sy^3 - exp(sx) + 3 == 0;
sol_vpa = vpasolve([eqn1, eqn2], [sx, sy], [p0x, p0y]);

x_vpa = double(sol_vpa.sx);
y_vpa = double(sol_vpa.sy);

fprintf('результат: x=%.10f, y=%.10f\n', x_vpa, y_vpa);
fprintf('проверка: f1=%.2e, f2=%.2e\n\n', ...
    0.1*y_vpa^3 - 5*x_vpa - y_vpa, ...
    y_vpa^3 - exp(x_vpa) + 3);

%метод простых итераций
fprintf('метод простых итераций:\n');

tau_val = 0.8;
%приведение к виду x = phi1(x,y), y = phi2(x,y)
psi1 = @(x,y) (0.1*y^3 - y) / 5;
psi2 = @(x,y) sign(exp(x) - 3) * abs(exp(x) - 3)^(1/3);

x_it = p0x; y_it = p0y;
max_it = 200; limit_tol = 1e-10;
it_count = 0; conv_ok = true;

fprintf('итерация |        x        |        y        |   ||j||  |   невязка\n');
fprintf('---------|-----------------|-----------------|---------|------------\n');

for k = 1:max_it
     x_tmp = psi1(x_it, y_it);
     y_tmp = psi2(x_it, y_it);
    
    % релаксация
    x_next = (1 - tau_val) * x_it + tau_val * x_tmp;
    y_next = (1 - tau_val) * y_it + tau_val * y_tmp;
    
    %вычисление якобиана для проверки сходимости
    eps_j = 1e-6;
    Jac11 = (psi1(x_it+eps_j, y_it) - psi1(x_it, y_it)) / eps_j;
    Jac12 = (psi1(x_it, y_it+eps_j) - psi1(x_it, y_it)) / eps_j;
    Jac21 = (psi2(x_it+eps_j, y_it) - psi2(x_it, y_it)) / eps_j;
    Jac22 = (psi2(x_it, y_it+eps_j) - psi2(x_it, y_it)) / eps_j;
    Norm_J = max(abs(eig([Jac11 Jac12; Jac21 Jac22])));
    
    %невязка
    fv1 = 0.1*y_it^3 - 5*x_it - y_it;
    fv2 = y_it^3 - exp(x_it) + 3;
    res_norm = norm([fv1; fv2]);
    
    fprintf('%9d | %.13f | %.13f | %7.4f | %.2e\n', ...
            k, x_next, y_next, Norm_J, res_norm);
    
    %проверка условия сходимости
    if Norm_J >= 1
        conv_ok = false;
    end
    
    %проверка остановки
    if max(abs([x_next - x_it, y_next - y_it])) < limit_tol && res_norm < limit_tol
        it_count = k;
        break;
    end
    
    x_it = x_next;
    y_it = y_next;
    it_count = k;
end

fprintf('\nрезультат: x=%.10f, y=%.10f\n', x_it, y_it);
fprintf('итераций: %d\n', it_count);
if conv_ok
    fprintf('условие сходимости выполнено на всех итерациях\n\n');
else
    fprintf('условие сходимости нарушено\n\n');
end

%метод зейделя
fprintf('метод зейделя:\n');

x_zei = p0x; y_zei = p0y;
max_it = 200; limit_tol = 1e-10;
it_zei = 0; conv_zei_ok = true;

fprintf('итерация |        x        |        y        |   ||j||  |   невязка\n');
fprintf('---------|-----------------|-----------------|---------|------------\n');

for k = 1:max_it
    %отличие зейделя: при вычислении y используем уже обновлённый x
    x_upd = psi1(x_zei, y_zei);
    y_upd = psi2(x_upd, y_zei);
    
    %вычисление якобиана для проверки сходимости
    eps_j = 1e-6;
    J11 = (psi1(x_zei+eps_j, y_zei) - psi1(x_zei, y_zei)) / eps_j;
    J12 = (psi1(x_zei, y_zei+eps_j) - psi1(x_zei, y_zei)) / eps_j;
    J21 = (psi2(x_zei+eps_j, y_zei) - psi2(x_zei, y_zei)) / eps_j;
    J22 = (psi2(x_zei, y_zei+eps_j) - psi2(x_zei, y_zei)) / eps_j;
    Norm_JZ = max(abs(eig([J11 J12; J21 J22])));
    
    %невязка
    rf1 = 0.1*y_zei^3 - 5*x_zei - y_zei;
    rf2 = y_zei^3 - exp(x_zei) + 3;
    res_zei = norm([rf1; rf2]);
    
    fprintf('%9d | %.13f | %.13f | %7.4f | %.2e\n', ...
            k, x_upd, y_upd, Norm_JZ, res_zei);
    
    %проверка условия сходимости
    if Norm_JZ >= 1
        conv_zei_ok = false;
    end
    
    %проверка остановки
    if max(abs([x_upd - x_zei, y_upd - y_zei])) < limit_tol && res_zei < limit_tol
        it_zei = k;
        break;
    end
    
    x_zei = x_upd;
    y_zei = y_upd;
    it_zei = k;
end

fprintf('\nрезультат метода зейделя: x = %.10f, y = %.10f, итераций: %d\n\n', x_zei, y_zei, it_zei);

%метод ньютона
fprintf('метод ньютона:\n');

fA = @(x,y) 0.1*y^3 - 5*x - y;
fB = @(x,y) y^3 - exp(x) + 3;

%аналитический якобиан
dfA_dx = @(x,y) -5;
dfA_dy = @(x,y) 0.3*y^2 - 1;
dfB_dx = @(x,y) -exp(x);
dfB_dy = @(x,y) 3*y^2;

x_nt = p0x; y_nt = p0y;
max_it_n = 100;
it_nt = 0;
lmbd = 1;

fprintf('итерация |        x        |        y        |   невязка  |    det(j)\n');
fprintf('---------|-----------------|-----------------|------------|------------\n');

for k = 1:max_it_n
    F_v = [fA(x_nt, y_nt); fB(x_nt, y_nt)];
    res_nt = norm(F_v);
    
    J_mat = [dfA_dx(x_nt, y_nt), dfA_dy(x_nt, y_nt);
             dfB_dx(x_nt, y_nt), dfB_dy(x_nt, y_nt)];
    det_J = det(J_mat);
    
    fprintf('%9d | %.13f | %.13f | %.2e | %.4e\n', ...
            k, x_nt, y_nt, res_nt, det_J);
    
    %проверка остановки
    if res_nt < 1e-12
        it_nt = k;
        break;
    end
    
    %проверка условия сходимости: якобиан не должен быть вырожден
    if abs(det_J) < 1e-8
        fprintf('якобиан близок к вырожденному, используется метод с релаксацией\n');
        dx_v = -J_mat \ F_v;
        lmbd = 0.5;
        x_nt = x_nt + lmbd * dx_v(1);
        y_nt = y_nt + lmbd * dx_v(2);
    else
        dx_v = -J_mat \ F_v;
        x_nt = x_nt + dx_v(1);
        y_nt = y_nt + dx_v(2);
    end
    
    it_nt = k;
end

fprintf('\nрезультат метода ньютона: x = %.10f, y = %.10f, итераций: %d\n\n', x_nt, y_nt, it_nt);

%финальный график
root_x = Res_fsolve(1);
root_y = Res_fsolve(2);

figure(2);
Lf1 = ezplot('0.1*y^3 - 5*x - y', [root_x-1, root_x+1, root_y-1, root_y+1]);
set(Lf1, 'Color', 'b', 'LineWidth', 2);
hold on;
Lf2 = ezplot('y^3 - exp(x) + 3', [root_x-1, root_x+1, root_y-1, root_y+1]);
set(Lf2, 'Color', 'r', 'LineWidth', 2);

%начальная точка - синяя точка
plot(p0x, p0y, 'bo', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'b');
%корень - красная звезда
plot(root_x, root_y, 'r*', 'MarkerSize', 15, 'LineWidth', 2);

xlabel('x'); ylabel('y');
title('решение системы нелинейных уравнений');
legend('f_1 = 0', 'f_2 = 0', 'начальная точка (синяя)', 'корень (красная *)', 'Location', 'best');
grid on; axis equal;

text(root_x+0.05, root_y+0.05, sprintf('(%.4f, %.4f)', root_x, root_y), ...
     'Color', 'red', 'FontWeight', 'bold');
 
%проверка подстановкой
fprintf('проверка:\n');
fprintf('подстановка корня (%.10f, %.10f) в исходные уравнения:\n', root_x, root_y);
fprintf('0.1*y^3 - 5*x - y = %.2e\n', 0.1*root_y^3 - 5*root_x - root_y);
fprintf('y^3 - exp(x) + 3 = %.2e\n\n', root_y^3 - exp(root_x) + 3);

%сводная таблица
fprintf('\nсводная таблица:\n');
fprintf('%-25s | %-12s | %-12s | %-10s\n', 'метод', 'x', 'y', 'итерации');
fprintf('--------------------------|--------------|--------------|------------\n');
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'fsolve', Res_fsolve(1), Res_fsolve(2), out_info.iterations);
fprintf('%-25s | %-12.8f | %-12.8f | %-10s\n', 'vpasolve', x_vpa, y_vpa, '-');
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'простые итерации', x_it, y_it, it_count);
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'метод зейделя', x_zei, y_zei, it_zei);
fprintf('%-25s | %-12.8f | %-12.8f | %-10d\n', 'метод ньютона', x_nt, y_nt, it_nt);