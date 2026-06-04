function App_Manual % Начало главной функции ручного расчета
    % Окно программы (полная копия интерфейса первой программы)
    hFig = figure('Name', 'Интерполяция и Аппроксимация (РУЧНАЯ МАТЕМАТИКА)', ...
        'Units', 'normalized', 'Position', [0.15 0.1 0.8 0.8], 'Color', [0.95 0.95 0.95], ...
        'MenuBar', 'none', 'NumberTitle', 'off');

    % Создание панели настроек
    panelControl = uipanel('Parent', hFig, 'Title', 'РУЧНОЙ РАСЧЕТ', ...
        'Units', 'normalized', 'Position', [0.01 0.02 0.28 0.96], 'FontWeight', 'bold');
    set(panelControl, 'Units', 'pixels'); % Фиксируем пиксели
    pPos = get(panelControl, 'Position'); panelWidth = pPos(3);

    % Ввод данных (те же дефолты: sin(x)+0.5x и веса)
    uicontrol(panelControl, 'Style', 'text', 'String', '1. Функция f(x):', ...
        'Position', [10 740 panelWidth-25 20], 'HorizontalAlignment', 'left');
    editFunction = uicontrol(panelControl, 'Style', 'edit', 'String', 'sin(x) + 0.5*x', ...
        'Position', [10 715 panelWidth-25 25], 'BackgroundColor', 'white');

    uicontrol(panelControl, 'Style', 'text', 'String', '2. Узлы X:', ...
        'Position', [10 685 panelWidth-25 20], 'HorizontalAlignment', 'left');
    editXNodes = uicontrol(panelControl, 'Style', 'edit', 'String', '0, 1.5, 3, 4.5, 6, 8, 10', ...
        'Position', [10 660 panelWidth-25 25], 'BackgroundColor', 'white');

    uicontrol(panelControl, 'Style', 'text', 'String', '3. Узлы Y (пусто = по f(x)):', ...
        'Position', [10 630 panelWidth-25 20], 'HorizontalAlignment', 'left');
    editYNodes = uicontrol(panelControl, 'Style', 'edit', 'String', '', ...
        'Position', [10 605 panelWidth-25 25], 'BackgroundColor', 'white');

    uicontrol(panelControl, 'Style', 'text', 'String', '4. Веса (точка 4 = приоритет):', ...
        'Position', [10 575 panelWidth-25 20], 'HorizontalAlignment', 'left', 'ForegroundColor', [0.7 0 0]);
    editWeights = uicontrol(panelControl, 'Style', 'edit', 'String', '1, 1, 1, 25, 1, 1, 1', ...
        'Position', [10 550 panelWidth-25 25], 'BackgroundColor', 'white');

    % Блок интерполяции Manual
    panelI = uipanel(panelControl, 'Title', 'Интерполяция', 'Units', 'pixels', ...
        'Position', [5 315 panelWidth-15 175], 'BackgroundColor', [1 0.9 0.9]);
    menuInterp = uicontrol(panelI, 'Style', 'popupmenu', ...
        'String', {'Линейная', 'Сплайн (СЛАУ)', 'Ньютон (Разности)', 'Лагранж (Циклы)'}, ...
        'Position', [10 115 panelWidth-40 25]);
    uicontrol(panelI, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ...
        'Position', [10 65 panelWidth-40 40], 'BackgroundColor', [0.8 0.2 0.2], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', 'Callback', @(~,~) runManual('interp'));
    uicontrol(panelI, 'Style', 'pushbutton', 'String', 'Матрица разностей', ...
        'Position', [10 15 panelWidth-40 40], 'Callback', @(~,~) callbackShowDiff());

    % Блок аппроксимации Manual
    panelA = uipanel(panelControl, 'Title', 'Аппроксимация', 'Units', 'pixels', ...
        'Position', [5 125 panelWidth-15 185], 'BackgroundColor', [1 1 0.8]);
    menuApprox = uicontrol(panelA, 'Style', 'popupmenu', ...
        'String', {'МНК Полином', 'Взвешенный МНК', 'Экспоненциальная', 'Степенная', 'Логарифмическая'}, ...
        'Position', [10 140 panelWidth-40 25]);
    uicontrol(panelA, 'Style', 'text', 'String', 'Степень m:', 'Position', [10 110 80 20], 'HorizontalAlignment', 'left');
    editDegree = uicontrol(panelA, 'Style', 'edit', 'String', '2', 'Position', [95 110 50 25], 'BackgroundColor', 'white');
    uicontrol(panelA, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ...
        'Position', [10 55 panelWidth-40 45], 'BackgroundColor', [0.5 0.5 0], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', 'Callback', @(~,~) runManual('approx'));
    uicontrol(panelA, 'Style', 'pushbutton', 'String', 'Найти оптимальное m', ...
        'Position', [10 10 panelWidth-40 35], 'Callback', @(~,~) callbackBestM());

    % Окно вывода ошибок
    textResultBox = uicontrol(panelControl, 'Style', 'edit', 'Max', 2, 'Min', 0, ...
        'Position', [10 10 panelWidth-25 110], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');

    set(panelControl, 'Units', 'normalized');
    axesMain = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.33 0.1 0.64 0.85]); grid on;

    % --- СБОР ДАННЫХ ---
    function [nodesX, nodesY, weights] = getDataFromUI()
        nodesX = str2num(get(editXNodes, 'String')); % Получаем X
        yIn = get(editYNodes, 'String'); % Получаем Y
        if isempty(yIn) % Если пусто, считаем по формуле
            f_a = eval(['@(x)' get(editFunction, 'String')]); nodesY = f_a(nodesX);
        else
            nodesY = str2num(yIn);
        end
        wIn = get(editWeights, 'String'); % Получаем веса
        if isempty(wIn), weights = ones(size(nodesX)); else, weights = str2num(wIn); end
    end

    % --- ГЛАВНЫЙ РАСЧЕТНЫЙ ЦИКЛ (РУЧНЫЕ АЛГОРИТМЫ) ---
    function runManual(mode)
        try
            [nX, nY, nW] = getDataFromUI(); % Считываем данные
            xf = linspace(min(nX), max(nX), 500); % Сетка графика
            numNodes = length(nX); % Число узлов
            yf = zeros(size(xf)); % Будущая кривая
            
            if strcmp(mode, 'interp') % Раздел интерполяции
                vI = get(menuInterp, 'Value'); % Метод
                
                if vI == 4 % ЛАГРАНЖ (РУКАМИ)
                    for k = 1:500 % Каждая точка графика
                        sumL = 0; % Накопление суммы
                        for i = 1:numNodes % Сумма по узлам
                            polyL = 1; % Начальное произведение
                            for j = 1:numNodes % Внутренний цикл
                                if i ~= j % Пропуск i-го индекса
                                    polyL = polyL * (xf(k) - nX(j)) / (nX(i) - nX(j)); % Перемножаем скобки
                                end
                            end
                            sumL = sumL + nY(i) * polyL; % Итоговое значение для i-го узла
                        end
                        yf(k) = sumL; % Записываем в итоговый массив
                    end
                elseif vI == 3 % НЬЮТОН (РАЗДЕЛЕННЫЕ РАЗНОСТИ РУКАМИ)
                    diffsTable = nY(:); % Первый столбец - значения Y
                    for j = 2:numNodes % Идем по порядку разности
                        for i = numNodes:-1:j % Считаем элементы снизу вверх
                            diffsTable(i) = (diffsTable(i) - diffsTable(i-1)) / (nX(i) - nX(i-j+1)); % Формула
                        end
                    end
                    for k = 1:500 % Считаем значение полинома Ньютона
                        val = diffsTable(numNodes); % Начинаем с последнего коэффициента
                        for i = numNodes-1:-1:1 % Двигаемся к первому по схеме Горнера
                            val = val * (xf(k) - nX(i)) + diffsTable(i); % Накопление результата
                        end
                        yf(k) = val; % Сохраняем в массив графика
                    end
                elseif vI == 1 % ЛИНЕЙНАЯ ИНТЕРПОЛЯЦИЯ РУКАМИ
                    for k = 1:500
                        idx = find(nX <= xf(k), 1, 'last'); % Ищем интервал для точки xf
                        if isempty(idx), idx = 1; end % Края
                        if idx < numNodes
                            yf(k) = nY(idx) + (nY(idx+1)-nY(idx)) * (xf(k)-nX(idx)) / (nX(idx+1)-nX(idx)); % Линейная формула
                        else, yf(k) = nY(numNodes); end
                    end
                elseif vI == 2 % КУБИЧЕСКИЙ СПЛАЙН (РЕШЕНИЕ СЛАУ РУКАМИ)
                    h_i = diff(nX); % Расстояния h
                    slopes = diff(nY)./h_i; % Наклоны сегментов
                    Amat = zeros(numNodes); % Матрица системы уравнений
                    Bvec = zeros(numNodes,1); % Вектор свободных членов
                    for i = 2:numNodes-1 % Заполняем СЛАУ для внутренних моментов
                        Amat(i, i-1) = h_i(i-1); % Значение h слева
                        Amat(i, i) = 2*(h_i(i-1) + h_i(i)); % Центральный коэффициент
                        Amat(i, i+1) = h_i(i); % Значение h справа
                        Bvec(i) = 3 * (slopes(i) - slopes(i-1)); % Свободный член (условие гладкости)
                    end
                    Amat(1,1)=1; Amat(numNodes, numNodes)=1; % Граничные условия
                    moments = Amat \ Bvec; % Решаем СЛАУ
                    for k = 1:500 % Собираем сплайн из сегментов
                        id = find(nX <= xf(k), 1, 'last'); if isempty(id) || id >= numNodes, id = numNodes-1; end
                        dx = xf(k) - nX(id); % Смещение x-xi
                        yf(k) = nY(id) + slopes(id)*dx + moments(id)*dx^2 + (moments(id+1)-moments(id))/(3*h_i(id))*dx^3; % Многочлен
                    end
                end
            else % Раздел аппроксимации
                m_a = str2double(get(editDegree, 'String')); % Степень m
                idxA = get(menuApprox, 'Value'); % Номер метода
                
                % ФИКС НУЛЕЙ (Адаптивный сдвиг для логарифмов и экспонент)
                Xs = nX; Ys = nY; xfc = xf;
                if idxA >= 3
                    if any(nX <= 0), Xs = nX + 0.1; xfc = xf + 0.1; end
                    if any(nY <= 0), Ys = nY + 0.1; end
                end

                if idxA <= 2 % РУЧНОЙ МНК
                    calcW = ones(1, numNodes); if idxA == 2, calcW = nW; end % Берем веса
                    V = zeros(numNodes, m_a+1); % Матрица Вандермонда руками
                    for r = 1:numNodes % Строки узлов
                        for c = 0:m_a % Столбцы степеней
                            V(r, c+1) = nX(r)^c; % Наполнение: x^0, x^1...
                        end
                    end
                    Wmat = diag(calcW); % Диагональная матрица из весов
                    coeffs = (V' * Wmat * V) \ (V' * Wmat * nY'); % РЕШЕНИЕ (V'WV)c = V'WY
                    for k = 1:500 % Считаем значения для графика по формуле полинома
                        ps = 0; for c = 0:m_a, ps = ps + coeffs(c+1) * xf(k)^c; end
                        yf(k) = ps;
                    end
                elseif idxA == 3 % РУЧНАЯ ЭКСПОНЕНТА
                    VM = [ones(numNodes,1), Xs(:)]; % Матрица для первой степени
                    C = (VM' * VM) \ (VM' * log(abs(Ys(:)))); % Решаем МНК для логарифма Y
                    yf = exp(C(1)) * exp(C(2) * xfc); % Обратный перевод в экспоненту
                elseif idxA == 4 % РУЧНАЯ СТЕПЕННАЯ
                    VM = [ones(numNodes,1), log(abs(Xs(:)))]; % Линейная матрица логарифмов
                    C = (VM' * VM) \ (VM' * log(abs(Ys(:)))); % Решаем МНК
                    yf = exp(C(1)) * xfc.^C(2); % Обратно в степень
                elseif idxA == 5 % РУЧНАЯ ЛОГАРИФМИЧЕСКАЯ
                    VM = [ones(numNodes,1), log(abs(Xs(:)))]; % Матрица
                    C = (VM' * VM) \ (VM' * Ys(:)); % Решаем МНК: Y = a + b*ln(x)
                    yf = C(1) + C(2) * log(abs(xfc)); % Результат
                end
            end
            
            % Финальная отрисовка
            cla(axesMain); plot(axesMain, xf, yf, 'g', 'LineWidth', 2); hold on;
            scatter(axesMain, nX, nY, 80, 'k', 'filled'); grid on;
            ya = zeros(1,numNodes); for i=1:numNodes, [~,minI]=min(abs(xf-nX(i))); ya(i)=yf(minI); end
            rmse = sqrt(mean((nY - ya).^2)); % Погрешность
            set(textResultBox, 'String', sprintf('РУЧНОЙ РАСЧЕТ\nRMSE: %.4e', rmse));
        catch ME, errordlg(ME.message); end
    end

    % Функции таблицы и подбора m абсолютно идентичныBuilt-in версии, но написаны циклом
    function callbackShowDiff(~, ~)
        [~, ny, ~] = getDataFromUI(); n = length(ny); T = zeros(n, n); T(:,1) = ny(:);
        for j = 2:n, for i = 1:n-j+1, T(i, j) = T(i+1, j-1) - T(i, j-1); end; end
        f = figure('Name','Разности'); uitable(f, 'Data', T, 'Units', 'normalized', 'Position', [0 0 1 1]);
    end

    function callbackBestM(~, ~)
        [nx, ny, ~] = getDataFromUI(); bestM = 1; minErr = inf;
        for m_t = 1:length(nx)-1
            V = zeros(length(nx), m_t+1); for i=1:length(nx), for j=0:m_t, V(i, j+1)=nx(i)^j; end; end
            C = (V' * V) \ (V' * ny'); ya = (V * C)'; err = sqrt(mean((ny-ya).^2));
            if err < minErr, minErr = err; bestM = m_t; end
        end
        set(editDegree, 'String', num2str(bestM));
    end
end