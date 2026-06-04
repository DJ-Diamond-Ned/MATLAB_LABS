function App_BuiltIn
    % -- % -- % == ОКНО ПРОГРАММЫ == % -- % -- %
    hFig = figure('Name', 'Интерполяция и Аппроксимация (Встроенные методы)', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85], 'Color',...
        [0.95 0.95 0.95], 'MenuBar', 'none', 'NumberTitle', 'off');

    % % % % % -- % -- % == ПАНЕЛЬ УПРАВЛЕНИЯ В ЛЕВОЙ ЧАСТИ ОКНА == % -- % -- % % % % %
    panelControl = uipanel('Parent', hFig,'Title', 'НАСТРОЙКИ РАСЧЕТА', 'Units', 'normalized','Position', [0.01 0.02 0.28 0.96], 'FontWeight', 'bold', 'FontSize', 10);

    set(panelControl, 'Units', 'pixels');
    panelPos = get(panelControl, 'Position');
    pW = panelPos(3);

    % -- % -- % == ПОД-СЕКЦИЯ 1: <<< ВВОД ФУНКЦИИ ПОЛЬЗОВАТЕЛЯ >>> == % -- % -- %
    uicontrol(panelControl, 'Style', 'text', 'String', '1. Функция f(x):', 'Position', [10 740 pW-25 20], 'HorizontalAlignment', 'left');
    hFun = uicontrol(panelControl, 'Style', 'edit', 'String', 'sin(x) + 0.5*x', 'Position', [10 715 pW-25 25], 'BackgroundColor', 'white');

    % -- % -- % == ПОД-СЕКЦИЯ 2: <<< ВВОД УЗЛОВЫХ ЗНАЧЕНИЙ ИКС >>> == % -- % -- %
    uicontrol(panelControl, 'Style', 'text', 'String', '2. Узлы X (через запятую):', 'Position', [10 685 pW-25 20], 'HorizontalAlignment', 'left');
    hX = uicontrol(panelControl, 'Style', 'edit', 'String', '0, 1.5, 3, 4.5, 6, 8, 10', 'Position', [10 660 pW-25 25], 'BackgroundColor', 'white');

    % -- % -- % == ПОД-СЕКЦИЯ 3: <<< ВВОД УЗЛОВЫХ ЗНАЧЕНИЙ ИГРЕК >>> == % -- % -- %
    uicontrol(panelControl, 'Style', 'text', 'String', '3. Узлы Y (пусто = по функции):', ... 
        'Position', [10 630 pW-25 20], 'HorizontalAlignment', 'left'); 
    hY = uicontrol(panelControl, 'Style', 'edit', 'String', '', ...
        'Position', [10 605 pW-25 25], 'BackgroundColor', 'white');

    % -- % -- % == ПОД-СЕКЦИЯ 4: <<< ВВОД ВЕСОВ УЗЛОВЫХ ТОЧЕК >>> == % -- % -- %
    uicontrol(panelControl, 'Style', 'text', 'String', '4. Веса:', 'Position', [10 575 pW-25 20], 'HorizontalAlignment', 'left', 'ForegroundColor', [0 0 0]);
    hW = uicontrol(panelControl, 'Style', 'edit', 'String', '1, 1, 1, 25, 1, 1, 1', 'Position', [10 550 pW-25 25], 'BackgroundColor', 'white');

    % % % -- % -- % == СИНЯЯ СЕКЦИЯ ИНТЕРПОЛЯЦИИ == % -- % -- % % %
    pInterp = uipanel(panelControl, 'Title', 'ИНТЕРПОЛЯЦИЯ', 'Units', 'pixels', 'Position', [5 365 pW-15 175], 'BackgroundColor', [0.9 0.95 1]);
    % -- % -- % == ДРОП-БОКС С МЕТОДАМИ ИНТЕРПОЛЯЦИИ == % -- % -- %
    hMethI = uicontrol(pInterp, 'Style', 'popupmenu', 'String', {'Linear (Линейная)', 'Spline (Сплайн)', 'PCHIP (Эрмитова)', 'Lagrange (Лагранж)'}, ...
        'Position', [10 125 pW-40 25]);
    % -- % -- % == КНОПКА ЗАПУСКА ИНТЕРПОЛЯЦИИ == % -- % -- %
    uicontrol(pInterp, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ...
        'Position', [10 75 pW-40 40], 'BackgroundColor', [0.2 0.4 0.8], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', 'Callback', @(~,~) runProcess('interp'));
    % -- % -- % == КНОПКА ВЫВОДА ТАБЛИЦЫ КОНЕЧНЫХ РАЗНОСТЕЙ == % -- % -- %
    uicontrol(pInterp, 'Style', 'pushbutton', 'String', 'Таблица разностей', ...
        'Position', [10 25 pW-40 40], 'Callback', @(~,~) showDiffWindow());

    % % % -- % -- % == ЗЕЛЁНАЯ СЕКЦИЯ АППРОКСИМАЦИИ == % -- % -- % % %
    pApprox = uipanel(panelControl, 'Title', 'АППРОКСИМАЦИЯ', 'Units', 'pixels', ...
        'Position', [5 170 pW-15 190], 'BackgroundColor', [0.9 1 0.9]);
    % -- % -- % == ДРОП-БОКС С МЕТОДАМИ АППРОКСИМАЦИИ == % -- % -- %
    hMethA = uicontrol(pApprox, 'Style', 'popupmenu', ...
        'String', {'МНК Полином', 'Взвешенный МНК', 'Экспонента', 'Степенная', 'Логарифм'}, ...
        'Position', [10 145 pW-40 25]);
    % -- % -- % == ТЕКСТ-БОКС ПОЛЕ ДЛЯ ВВОДА m (СТЕПЕНИ ПОЛИНОМА) == % -- % -- %
    uicontrol(pApprox, 'Style', 'text', 'String', 'Степень m:', 'Position', [10 115 80 20], 'HorizontalAlignment', 'left');
    hDeg = uicontrol(pApprox, 'Style', 'edit', 'String', '2', 'Position', [95 110 50 25], 'BackgroundColor', 'white');
    % -- % -- % == КНОПКА ЗАПУСКА АППРОКСИМАЦИИ == % -- % -- %
    uicontrol(pApprox, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ...
        'Position', [10 60 pW-40 45], 'BackgroundColor', [0.1 0.6 0.1], 'ForegroundColor', 'white', ...
        'FontWeight', 'bold', 'Callback', @(~,~) runProcess('approx'));
    % -- % -- % == КНОПКА ЗАПУСКА ПОИСКА ЛУЧШЕЙ СТЕПЕНИ ПОЛИНОМА == % -- % -- %
    uicontrol(pApprox, 'Style', 'pushbutton', 'String', 'Найти лучшую степень m', ...
        'Position', [10 15 pW-40 35], 'Callback', @(~,~) findBestDegree());

    % -- % -- % == НИЖНЯЯ ЧАСТЬ ДЛЯ ВЫВОДА RMSE И РЕЗУЛЬТАТОВ == % -- % -- %
    hResText = uicontrol(panelControl, 'Style', 'edit', 'Max', 2, 'Min', 0, ...
        'Position', [10 10 pW-25 150], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
        
    set(panelControl, 'Units', 'normalized');

    % % % -- % -- % == ОБЛАСТЬ ДЛЯ ОТРИСОВКИ ГРАФИКОВ ВВЕДЁННОЙ ПОЛЬЗОВАТЕЛЕМ ФУНКЦИИ
    hAxes = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.33 0.1 0.64 0.85]);
    grid(hAxes, 'on');

    % < > % ===== % -- ВНУТРЕННЯЯ ФУНКЦИЯ СБОРА ДАННЫХ (getDataFromUI) -- % ===== % < > %
    function [nodesX, nodesY, weights] = getDataFromUI()
        rawX = str2num(get(hX, 'String')); % Конвертируем текст из поля X в массив чисел
        nodesX = rawX(:); % Принудительно превращаем массив X в вертикальный столбец
        
        yInputText = get(hY, 'String'); % Получаем текст из поля ввода узлов Y
        if isempty(yInputText) % Если пользователь не ввел Y вручную
            funcString = get(hFun, 'String'); % Берем формулу функции f(x) из поля ввода
            fHandle = eval(['@(x)' funcString]); % Превращаем строку в исполняемую анонимную функцию MATLAB
            rawY = fHandle(nodesX); % Вычисляем значения Y для всех заданных точек X
            
        else % Если пользователь ввел значения Y через запятую
            rawY = str2num(yInputText); % Конвертируем строку Y в числовой массив
        end
        
        nodesY = rawY(:); % Принудительно превращаем массив Y в вертикальный столбец
        
        wInputText = get(hW, 'String'); % Считываем текст из поля весов для МНК
        if isempty(wInputText) % Если поле весов пустое
            weights = ones(size(nodesX)); % Создаем массив единиц такой же длины, как X
            
        else % Если веса введены пользователем
            rawW = str2num(wInputText); % Превращаем строку весов в числовой массив
            weights = rawW(:); % Принудительно превращаем массив весов в вертикальный столбец
        end
    end

    % < > % ===== % -- ГЛАВНАЯ ФУНКЦИЯ РАСЧЕТА (runProcess) -- % ===== % < > %
    function runProcess(mode) % Эта функция вызывается кнопками РАССЧИТАТЬ
        try % Начало блока обработки ошибок (чтобы программа не вылетала при неверном вводе)
            [nodesX, nodesY, weights] = getDataFromUI(); % Вызываем функцию получения актуальных данных из полей
            fineX = linspace(min(nodesX), max(nodesX), 500)'; % Генерируем 500 точек для рисования (столбец)
            
            if strcmp(mode, 'interp') % --- ЕСЛИ НАЖАТА КНОПКА ИНТЕРПОЛЯЦИИ ---
                methodIdx = get(hMethI, 'Value'); % Получаем номер выбранного метода из списка
                
                if methodIdx == 1 % Выбрана Линейная
                    fineY = interp1(nodesX, nodesY, fineX, 'linear'); % Встроенная функция линейной интерполяции
                    label = 'Linear'; % Запоминаем название для вывода
                    
                elseif methodIdx == 2 % Выбран Сплайн
                    fineY = interp1(nodesX, nodesY, fineX, 'spline'); % Встроенная функция кубического сплайна
                    label = 'Spline'; % Запоминаем название
                    
                elseif methodIdx == 3 % Выбрана Эрмитова (PCHIP)
                    fineY = interp1(nodesX, nodesY, fineX, 'pchip'); % Встроенная функция Эрмитова сплайна
                    label = 'PCHIP'; % Запоминаем название
                    
                elseif methodIdx == 4 % Выбран Лагранж
                    polyCoeffs = polyfit(nodesX, nodesY, length(nodesX)-1); % Полином степени N-1 (Лагранж)
                    fineY = polyval(polyCoeffs, fineX); % Считаем значения полинома на сетке fineX
                    label = 'Lagrange'; % Запоминаем название
                end
                
            else % --- ЕСЛИ НАЖАТА КНОПКА АППРОКСИМАЦИИ ---
                mDegree = str2double(get(hDeg, 'String')); % Считываем введенную степень полинома m
                approxIdx = get(hMethA, 'Value'); % Получаем индекс выбранного метода аппроксимации
                shiftVal = 0.05; % Маленькое смещение для защиты от логарифма нуля
                
                if approxIdx == 1 % МНК Полином
                    coeffsPoly = polyfit(nodesX, nodesY, mDegree); % Считаем коэффициенты полинома степени m
                    fineY = polyval(coeffsPoly, fineX); % Считаем значения линии графика
                    label = 'МНК Полином'; % Название метода
                    
                elseif approxIdx == 2 % Взвешенный МНК
                    % Формируем матрицу Вандермонда вручную для точности
                    VMat = zeros(length(nodesX), mDegree + 1); % Пустая матрица размером N на m+1
                    for p = 0:mDegree % Цикл по степеням от 0 до m
                        VMat(:, p + 1) = nodesX.^p; % Наполняем столбцы матрицы: x^0, x^1, x^2...
                    end
                    
                    coeffsW = lscov(VMat, nodesY, weights); % Решаем задачу WLS через встроенную lscov
                    fineY = zeros(size(fineX)); % Пустой массив для отрисовки
                    
                    for p = 0:mDegree % Цикл сборки результирующего полинома по найденным коэффициентам
                        fineY = fineY + coeffsW(p + 1) * fineX.^p; % Накопление суммы c_i * x^i
                    end 
                    
                    label = 'Взвешенный МНК'; % Название метода
                    
                elseif approxIdx == 3 % Экспонента y = a * e^(bx)
                    % Линеаризация: решаем ln(y+shift) = ln(a) + bx через polyfit 1-й степени
                    pExp = polyfit(nodesX, log(abs(nodesY) + shiftVal), 1); % Коэффициенты логарифма
                    fineY = exp(pExp(2)) * exp(pExp(1) * fineX) - shiftVal; % Обратное преобразование в формулу
                    label = 'Экспонента'; % Название
                    
                elseif approxIdx == 4 % Степенная y = a * x^b
                    % Линеаризация: ln(y+shift) по ln(x+shift)
                    pPow = polyfit(log(nodesX + shiftVal), log(abs(nodesY) + shiftVal), 1); % Считаем коэффициенты
                    fineY = exp(pPow(2)) * (fineX + shiftVal).^pPow(1) - shiftVal; % Считаем результат
                    label = 'Степенная'; % Название
                    
                elseif approxIdx == 5 % Логарифм y = a * ln(x+shift) + b
                    % Линеаризация: Y по ln(x+shift)
                    pLog = polyfit(log(nodesX + shiftVal), nodesY, 1); % Считаем коэффициенты МНК
                    fineY = pLog(1) * log(fineX + shiftVal) + pLog(2); % Считаем результат
                    label = 'Логарифм'; % Название
                end
            end
            
            
            % --- ВИЗУАЛИЗАЦИЯ И ВЫВОД РЕЗУЛЬТАТОВ ---
            cla(hAxes); % Полностью очищаем область графика перед рисованием новой кривой
            plot(hAxes, fineX, fineY, 'r', 'LineWidth', 2); % Рисуем красную расчетную кривую
            hold(hAxes, 'on'); % Включаем режим наложения (hold on), чтобы добавить точки
            scatter(hAxes, nodesX, nodesY, 80, 'b', 'filled'); % Рисуем исходные узлы синими точками
            grid on; % Показываем сетку координат
            title(hAxes, ['Результат: ' label]); % Ставим заголовок с названием метода
            
            % Считаем отклонение (RMSE) в точках узлов для оценки точности
            yNodesCheck = interp1(fineX, fineY, nodesX, 'linear', 'extrap'); % Получаем значения с кривой в узлах
            rmseResult = sqrt(mean((nodesY - yNodesCheck).^2)); % Формула среднеквадратичной ошибки (RMSE)
            maxErrorVal = max(abs(nodesY - yNodesCheck)); % Формула максимального отклонения в узле
            
            % Формируем текст для окна вывода результатов
            resString = sprintf('МЕТОД: %s\nRMSE (Погрешность): %.4e\nMAX ERROR (Макс. откл): %.4e', ...
                label, rmseResult, maxErrorVal); % Создаем строку
            set(hResText, 'String', resString); % Выводим текст в интерфейс программы
            
        catch ME
            errordlg(['Ошибка в расчете: ' ME.message]); % Показываем всплывающее окно с текстом ошибки
        end
    end

    % --- ТАБЛИЦА КОНЕЧНЫХ РАЗНОСТЕЙ ---
    function showDiffWindow() % Функция для кнопки "Таблица разностей"
        [~, currentNodesY, ~] = getDataFromUI(); % Получаем актуальные значения Y из интерфейса
        totalN = length(currentNodesY); % Считаем количество точек
        diffTable = zeros(totalN, totalN); % Создаем пустую квадратную матрицу для разностей
        diffTable(:, 1) = currentNodesY(:); % Записываем значения Y в первый столбец матрицы
        for col = 2:totalN % Цикл по порядку разности (2-я, 3-я и так далее)
            for row = 1:totalN - col + 1 % Цикл заполнения строк в текущем столбце
                % Конечная разность = значение ниже минус значение выше из предыдущего столбца
                diffTable(row, col) = diffTable(row+1, col-1) - diffTable(row, col-1); 
            end
        end
        hDiffFig = figure('Name', 'Таблица конечных разностей', 'MenuBar', 'none', 'NumberTitle', 'off', 'Position', [200 200 600 400]); % Создаем новое окно
        uitable('Parent', hDiffFig, 'Data', diffTable, 'Position', [20 20 560 360]);  % Отрисовываем таблицу
    end

    % --- АВТОПОДБОР СТЕПЕНИ m ---
    function findBestDegree() % Функция для поиска наилучшей степени m по RMSE
        [currX, currY, ~] = getDataFromUI(); % Считываем текущие узлы X и Y
        bestDegreeFound = 1; % По умолчанию считаем лучшей степень 1
        minimumRMSE = inf; % Изначально ставим ошибку равной бесконечности
        for m_test = 1:length(currX)-1 % Перебираем все возможные степени полинома от 1 до N-1
            p_test = polyfit(currX, currY, m_test); % Считаем МНК для текущей степени m_test
            y_test = polyval(p_test, currX); % Считаем значения полинома в точках узлов
            currentRMSE = sqrt(mean((currY - y_test).^2)); % Вычисляем RMSE отклонение
            if currentRMSE < minimumRMSE % Если найденная ошибка меньше самого лучшего результата
                minimumRMSE = currentRMSE; % Запоминаем новую минимальную ошибку
                bestDegreeFound = m_test; % Запоминаем текущую степень m как лучшую
            end
        end
        set(hDeg, 'String', num2str(bestDegreeFound)); % Обновляем значение степени m в интерфейсе
    end
end