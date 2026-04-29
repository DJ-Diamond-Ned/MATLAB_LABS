A = [2.34 -1.42 -0.54 0.21;
     1.44 -0.53 1.43 -1.27;
     0.63 -1.32 -0.65 1.43;
     0.54 0.88 -0.67 -2.38];

B = [0.66;-1.44;0.94;0.73];

% Характеристики
detA = det(A)
rankA = rank(A)
normA = norm(A)

% LU разложение
[L,U] = lu(A);

y = L\B;
x1 = U\y

% Через linsolve
x2 = linsolve(A,B)

% Проверка
check1 = A*x1 - B
check2 = A*x2 - B

% Обусловленность
condA = cond(A)