A = [6 -1 -1;
     1 -2 3;
     3 4 4];

B = [0;1;-1];

% Характеристики матрицы
detA = det(A)
rankA = rank(A)
normA = norm(A)

% Решение методом обратной матрицы
x1 = inv(A)*B

% Решение через linsolve
x2 = linsolve(A,B)

% Проверка точности
check1 = A*x1 - B
check2 = A*x2 - B

% Обусловленность
condA = cond(A) 