clear; clc; close all;

% параметры интегрирования
h = 0.1;
x_end1 = 1;
x_end2 = 3;
x1 = 0:h:x_end1;
x2 = 0:h:x_end2;
n1 = length(x1);
n2 = length(x2);

f1_neg = @(x, y1, y2) y1*exp(-x^2) + x*y2;
f1_pos = @(x, y1, y2) y1*exp(x^2) + x*y2;
f2 = @(x, y1, y2) 3*x - y1 + 2*y2;

% априорная стратегия выбора шага
% оцениваем максимальное собственное число матрицы якоби на интервале
fprintf('априорная стратегия выбора шага:\n');

% пробные точки для оценки
x_sample = 0:0.1:3;
max_lambda = 0;
for i = 1:length(x_sample)
    % матрица якоби для системы
    J = [exp(-x_sample(i)^2), x_sample(i); -1, 2];
    eig_vals = eig(J);
    max_lambda = max(max_lambda, max(abs(eig_vals)));
end
% условие устойчивости явного метода эйлера
h_max_stable = 2 / max_lambda;
h_safe = h_max_stable / 2;  % берём с запасом

fprintf('  максимальное |lambda| = %.2f\n', max_lambda);
fprintf('  максимальный устойчивый шаг: h_max = %.4f\n', h_max_stable);
fprintf('  выбранный априорный шаг: h = %.4f\n', h_safe);
fprintf('  (по заданию используем h = 0.1, явные методы могут быть неустойчивы)\n\n');

% массив для хранения числа жёсткости на каждом шаге
stiffness1 = zeros(1, n1);

% явный метод эйлера
y1_eul = zeros(1,n1); y2_eul = zeros(1,n1);
y1_eul(1) = 1; y2_eul(1) = 1;  % начальные условия

for i = 1:n1-1
    % формула явного эйлера: y_{n+1} = y_n + h * f(x_n, y_n)
    y1_eul(i+1) = y1_eul(i) + h*f1_neg(x1(i), y1_eul(i), y2_eul(i));
    y2_eul(i+1) = y2_eul(i) + h*f2(x1(i), y1_eul(i), y2_eul(i));
    
    % вычисление числа жёсткости на текущем шаге
    % число жёсткости = |lambda_max| / |lambda_min|
    J = [exp(-x1(i)^2), x1(i); -1, 2];
    eig_vals = eig(J);
    stiffness1(i) = max(abs(eig_vals)) / min(abs(eig_vals));
end
stiffness1(end) = stiffness1(end-1);  % для последней точки

fprintf('явный метод эйлера:\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_eul(1), y2_eul(1));
fprintf('  x=0.5: y1=%.6f, y2=%.6f\n', y1_eul(6), y2_eul(6));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_eul(end), y2_eul(end));
fprintf('  число жёсткости на шаге x=0.5: %.2f\n\n', stiffness1(6));

% модифицированный метод эйлера
y1_mod = zeros(1,n1); y2_mod = zeros(1,n1);
y1_mod(1) = 1; y2_mod(1) = 1;

for i = 1:n1-1
    % первый наклон - в начале шага
    k1y1 = f1_neg(x1(i), y1_mod(i), y2_mod(i));
    k1y2 = f2(x1(i), y1_mod(i), y2_mod(i));
    
    % предсказание по явному эйлеру
    y1_pred = y1_mod(i) + h*k1y1;
    y2_pred = y2_mod(i) + h*k1y2;
    
    % второй наклон - в конце шага (на предсказанных значениях)
    k2y1 = f1_neg(x1(i)+h, y1_pred, y2_pred);
    k2y2 = f2(x1(i)+h, y1_pred, y2_pred);
    
    % итоговое значение - среднее арифметическое наклонов
    y1_mod(i+1) = y1_mod(i) + h/2*(k1y1 + k2y1);
    y2_mod(i+1) = y2_mod(i) + h/2*(k1y2 + k2y2);
end

fprintf('модифицированный метод эйлера:\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_mod(1), y2_mod(1));
fprintf('  x=0.5: y1=%.6f, y2=%.6f\n', y1_mod(6), y2_mod(6));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n\n', y1_mod(end), y2_mod(end));

% метод рунге-кутта 4 порядка
y1_rk4 = zeros(1,n1); y2_rk4 = zeros(1,n1);
y1_rk4(1) = 1; y2_rk4(1) = 1;

for i = 1:n1-1
    % коэффициент k1 - наклон в начале
    k1 = h * [f1_neg(x1(i), y1_rk4(i), y2_rk4(i)); 
              f2(x1(i), y1_rk4(i), y2_rk4(i))];
    
    % коэффициент k2 - наклон в середине (с использованием k1)
    k2 = h * [f1_neg(x1(i)+h/2, y1_rk4(i)+k1(1)/2, y2_rk4(i)+k1(2)/2);
              f2(x1(i)+h/2, y1_rk4(i)+k1(1)/2, y2_rk4(i)+k1(2)/2)];
    
    % коэффициент k3 - наклон в середине (с использованием k2)
    k3 = h * [f1_neg(x1(i)+h/2, y1_rk4(i)+k2(1)/2, y2_rk4(i)+k2(2)/2);
              f2(x1(i)+h/2, y1_rk4(i)+k2(1)/2, y2_rk4(i)+k2(2)/2)];
    
    % коэффициент k4 - наклон в конце (с использованием k3)
    k4 = h * [f1_neg(x1(i)+h, y1_rk4(i)+k3(1), y2_rk4(i)+k3(2));
              f2(x1(i)+h, y1_rk4(i)+k3(1), y2_rk4(i)+k3(2))];
    
    % итоговое значение - взвешенное среднее четырёх наклонов
    y1_rk4(i+1) = y1_rk4(i) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
    y2_rk4(i+1) = y2_rk4(i) + (k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6;
end

fprintf('метод рунге-кутта 4 порядка:\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_rk4(1), y2_rk4(1));
fprintf('  x=0.5: y1=%.6f, y2=%.6f\n', y1_rk4(6), y2_rk4(6));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n\n', y1_rk4(end), y2_rk4(end));

% эталонное решение через встроенный решатель matlab ode45
% ode45 - адаптивный метод рунге-кутта, используется как эталон
[x_ode1, y_ode1] = ode45(@(x,y) [y(1)*exp(-x^2) + x*y(2); 3*x - y(1) + 2*y(2)], x1, [1;1]);
y1_ode1 = interp1(x_ode1, y_ode1(:,1), x1);
y2_ode1 = interp1(x_ode1, y_ode1(:,2), x1);

fprintf('решение стандартным оператором matlab (ode45):\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_ode1(1), y2_ode1(1));
fprintf('  x=0.5: y1=%.6f, y2=%.6f\n', y1_ode1(6), y2_ode1(6));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n\n', y1_ode1(end), y2_ode1(end));

% графики для части 1
% график y1(x)
figure('Name', 'часть 1 y1(x)');
plot(x1, y1_eul, 'b-o', x1, y1_mod, 'g-s', x1, y1_rk4, 'r-^', x1, y1_ode1, 'k--', 'LineWidth',1.5);
xlabel('x'); ylabel('y1(x)'); title('часть 1: сравнение методов y1(x)');
legend('явный эйлер', 'модиф. эйлер', 'рунге-кутта 4', 'ode45', 'Location', 'best');
grid on;

% график y2(x)
figure('Name', 'часть 1 y2(x)');
plot(x1, y2_eul, 'b-o', x1, y2_mod, 'g-s', x1, y2_rk4, 'r-^', x1, y2_ode1, 'k--', 'LineWidth',1.5);
xlabel('x'); ylabel('y2(x)'); title('часть 1: сравнение методов y2(x)');
legend('явный эйлер', 'модиф. эйлер', 'рунге-кутта 4', 'ode45', 'Location', 'best');
grid on;

% график числа жёсткости
figure('Name', 'часть 1 число жёсткости');
plot(x1, stiffness1, 'b-', 'LineWidth',1.5);
xlabel('x'); ylabel('число жёсткости'); title('число жёсткости на каждом шаге (часть 1)');
grid on;

stiffness2 = zeros(1, n2);

% ----- явный метод эйлера для части 2 -----
y1_eul2 = zeros(1,n2); y2_eul2 = zeros(1,n2);
y1_eul2(1) = 1; y2_eul2(1) = 1;

for i = 1:n2-1
    y1_eul2(i+1) = y1_eul2(i) + h*f1_pos(x2(i), y1_eul2(i), y2_eul2(i));
    y2_eul2(i+1) = y2_eul2(i) + h*f2(x2(i), y1_eul2(i), y2_eul2(i));
    
    % вычисление числа жёсткости (теперь с e^(x^2))
    J = [exp(x2(i)^2), x2(i); -1, 2];
    eig_vals = eig(J);
    stiffness2(i) = max(abs(eig_vals)) / min(abs(eig_vals));
end
stiffness2(end) = stiffness2(end-1);

fprintf('явный метод эйлера:\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_eul2(1), y2_eul2(1));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_eul2(11), y2_eul2(11));
fprintf('  x=2:   y1=%.6f, y2=%.6f\n', y1_eul2(21), y2_eul2(21));
fprintf('  x=3:   y1=%.6f, y2=%.6f\n', y1_eul2(end), y2_eul2(end));
fprintf('  число жёсткости на шаге x=1: %.2f\n\n', stiffness2(11));

% неявный метод эйлера с символьным решением
% символьные переменные для аналитического решения системы
syms y1n y2n y1_old y2_old h_sym x_sym

% записываем неявные уравнения:
% y1_new = y1_old + h * (y1_new * exp(x^2) + x * y2_new)
% y2_new = y2_old + h * (3*x - y1_new + 2*y2_new)
eq1_sym = y1n == y1_old + h_sym*(y1n*exp(x_sym^2) + x_sym*y2n);
eq2_sym = y2n == y2_old + h_sym*(3*x_sym - y1n + 2*y2n);

% решаем систему аналитически относительно y1n и y2n
sol_sym = solve([eq1_sym eq2_sym], [y1n y2n]);

% преобразуем символьные решения в числовые функции
y1_next = matlabFunction(sol_sym.y1n, 'Vars', [y1_old, y2_old, h_sym, x_sym]);
y2_next = matlabFunction(sol_sym.y2n, 'Vars', [y1_old, y2_old, h_sym, x_sym]);

fprintf('неявный метод эйлера (символьное решение):\n');
fprintf('  получены аналитические формулы для y1_new и y2_new\n\n');

% применяем неявный метод с шагом 0.1
y1_imp = zeros(1,n2); y2_imp = zeros(1,n2);
y1_imp(1) = 1; y2_imp(1) = 1;

for i = 1:n2-1
    % используем символьные формулы для расчёта следующего шага
    y1_imp(i+1) = y1_next(y1_imp(i), y2_imp(i), h, x2(i)+h);
    y2_imp(i+1) = y2_next(y1_imp(i), y2_imp(i), h, x2(i)+h);
end

fprintf('неявный метод эйлера (символьное решение):\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_imp(1), y2_imp(1));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_imp(11), y2_imp(11));
fprintf('  x=2:   y1=%.6f, y2=%.6f\n', y1_imp(21), y2_imp(21));
fprintf('  x=3:   y1=%.6f, y2=%.6f\n\n', y1_imp(end), y2_imp(end));

% эталонное решение ode45 для части 2
[x_ode2, y_ode2] = ode45(@(x,y) [y(1)*exp(x^2) + x*y(2); 3*x - y(1) + 2*y(2)], x2, [1;1]);
y1_ode2 = interp1(x_ode2, y_ode2(:,1), x2);
y2_ode2 = interp1(x_ode2, y_ode2(:,2), x2);

fprintf('решение стандартным оператором matlab (ode45):\n');
fprintf('  x=0:   y1=%.6f, y2=%.6f\n', y1_ode2(1), y2_ode2(1));
fprintf('  x=1:   y1=%.6f, y2=%.6f\n', y1_ode2(11), y2_ode2(11));
fprintf('  x=2:   y1=%.6f, y2=%.6f\n', y1_ode2(21), y2_ode2(21));
fprintf('  x=3:   y1=%.6f, y2=%.6f\n\n', y1_ode2(end), y2_ode2(end));

% графики для части 2
% график y1(x)
figure('Name', 'часть 2 y1(x)');
plot(x2, y1_eul2, 'b-o', x2, y1_imp, 'm-s', x2, y1_ode2, 'k--', 'LineWidth',1.5);
xlabel('x'); ylabel('y1(x)'); title('часть 2: явный vs неявный эйлер y1(x)');
legend('явный эйлер', 'неявный эйлер (символьн.)', 'ode45', 'Location', 'best');
grid on;

% график y2(x)
figure('Name', 'часть 2 y2(x)');
plot(x2, y2_eul2, 'b-o', x2, y2_imp, 'm-s', x2, y2_ode2, 'k--', 'LineWidth',1.5);
xlabel('x'); ylabel('y2(x)'); title('часть 2: явный vs неявный эйлер y2(x)');
legend('явный эйлер', 'неявный эйлер (символьн.)', 'ode45', 'Location', 'best');
grid on;

% график числа жёсткости для части 2
figure('Name', 'часть 2 число жёсткости');
plot(x2, stiffness2, 'b-', 'LineWidth',1.5);
xlabel('x'); ylabel('число жёсткости'); title('число жёсткости на каждом шаге (часть 2)');
grid on;