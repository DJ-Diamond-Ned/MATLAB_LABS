A = [9.1 5.6 7.8;
     3.8 5.1 2.8;
     4.1 5.7 1.2];

B = [9.8;6.7;5.8];

% Характеристики матрицы
detA = det(A)
rankA = rank(A)
normA = norm(A)

% Метод Гаусса
x1 = A\B

% Через linsolve
x2 = linsolve(A,B)

% Проверка
check1 = A*x1 - B
check2 = A*x2 - B

% Обусловленность
condA = cond(A)