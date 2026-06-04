clc; clear;

% 1. Исходные данные
substances = {'Na2CO3', 'HNO3', 'NaNO3', 'H2O', 'CO2', 'CaO', 'Ca(NO3)2'};
% Атомная матрица (строки: Na, C, O, H, N, Ca)
A = [2 0 1 0 0 0 0;  % Na
     1 0 0 0 1 0 0;  % C
     3 3 3 1 2 1 6;  % O
     0 1 0 2 0 0 0;  % H
     0 1 1 0 0 0 2;  % N
     0 0 0 0 0 1 1]; % Ca

% Находим ранг и опорные столбцы
[~, pivots] = rref(A); 
rankA = rank(A);
free_cols = setdiff(1:7, pivots);

fprintf('--- Общие характеристики ---\n');
fprintf('Ранг основной матрицы A = %d\n', rankA);
fprintf('Норма матрицы A = %.2f\n\n', norm(A));

all_reactions = [];

% Цикл по свободным переменным (находим 2 независимых набора реакций)
for i = 1:length(free_cols)
    fprintf('--- Анализ набора реакций №%d ---\n', i);
    
    current_free = free_cols(i);
    
    row_indices = [2 3 4 5 6]; 
    B_sub = A(row_indices, pivots); 
    f_sub = -A(row_indices, current_free);
    
    % Вычисления для подматриц по критериям (det, rank, norm, cond)
    detB = det(B_sub);
    fprintf('Определитель подматрицы: %.2f\n', detB);
    fprintf('Ранг подматрицы: %d\n', rank(B_sub));
    fprintf('Норма подматрицы: %.2f\n', norm(B_sub));
    fprintf('Обусловленность подматрицы: %.2f\n', cond(B_sub));
    
    % Решение СЛАУ методом linsolve (Критерий 0.1)
    x_part = linsolve(B_sub, f_sub);
    
    % Сборка полного вектора коэффициентов
    res = zeros(7, 1);
    res(pivots) = x_part;
    res(current_free) = 1;
    
    % Безопасное приведение к целым числам (исправление ошибки lcm)
    [num, den] = rat(res, 1e-5);
    den = round(abs(den)); % Принудительно делаем целыми для lcm
    L = den(1);
    for k = 2:length(den)
        if den(k) ~= 0
            L = lcm(L, den(k));
        end
    end
    res_int = round(res * L);
    
    % Проверка точности Ax=0 (Критерий 0.1)
    fprintf('Точность (Ax=0): норма ошибки = %e\n\n', norm(A * res_int));
    
    all_reactions = [all_reactions, res_int];
end

% 4. Автоматизированный вывод уравнений (Бонус 0.3)
fprintf('--- Итоговые химические уравнения ---\n');
for r = 1:size(all_reactions, 2)
    coeffs = all_reactions(:, r);
    left = {}; right = {};
    for c = 1:length(coeffs)
        val = coeffs(c);
        if val > 0
            left{end+1} = sprintf('%d%s', val, substances{c});
        elseif val < 0
            right{end+1} = sprintf('%d%s', abs(val), substances{c});
        end
    end
    fprintf('Реакция %d: %s = %s\n', r, strjoin(left, ' + '), strjoin(right, ' + '));
end