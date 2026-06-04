function App_BuiltIn % Объявление главной функции, в которой живет всё приложение
    % Создаем основное графическое окно программы
    hFig = figure('Name', 'Интерполяция и Аппроксимация (Built-in)', ... % Заголовок
        'Units', 'normalized', ... % Координаты в процентах от экрана
        'Position', [0.1 0.1 0.8 0.8], ... % Окно: отступ 10%, размер 80%
        'Color', [0.95 0.95 0.95], ... % Цвет фона (светло-серый)
        'MenuBar', 'none', ... % Отключаем стандартное меню MATLAB
        'NumberTitle', 'off'); % Убираем надпись "Figure 1"

    % --- СОЗДАНИЕ ЛЕВОЙ ПАНЕЛИ УПРАВЛЕНИЯ ---
    panelControl = uipanel('Parent', hFig, ... % Панель внутри главного окна
        'Title', 'НАСТРОЙКИ ДАННЫХ', ... % Заголовок на рамке
        'Units', 'normalized', ... % В процентах
        'Position', [0.01 0.02 0.28 0.96], ... % Слева, почти на всю высоту
        'FontWeight', 'bold', ... % Жирный шрифт заголовка
        'FontSize', 10); % Размер шрифта 10

    % Переходим на пиксели для точного расположения элементов внутри панели
    set(panelControl, 'Units', 'pixels'); % Теперь Position задается в пикселях
    panelPos = get(panelControl, 'Position'); % Получаем реальные размеры панели
    pW = panelPos(3); % Сохраняем ширину панели в переменную pW

    % 1. Секция ввода функции f(x)
    uicontrol(panelControl, 'Style', 'text', 'String', '1. Функция f(x):', ... % Метка
        'Position', [10 740 pW-25 20], 'HorizontalAlignment', 'left'); % Слева
    hEditFun = uicontrol(panelControl, 'Style', 'edit', 'String', 'sin(x) + 0.5*x', ... % Поле ввода
        'Position', [10 715 pW-25 25], 'BackgroundColor', 'white'); % Белый фон

    % 2. Секция ввода узлов X
    uicontrol(panelControl, 'Style', 'text', 'String', '2. Узлы X (через запятую):', ... % Метка
        'Position', [10 685 pW-25 20], 'HorizontalAlignment', 'left'); % Слева
    hEditX = uicontrol(panelControl, 'Style', 'edit', 'String', '0, 1.5, 3, 4.5, 6, 8, 10', ... % Узлы X
        'Position', [10 660 pW-25 25], 'BackgroundColor', 'white'); % Белый фон

    % 3. Секция ввода узлов Y
    uicontrol(panelControl, 'Style', 'text', 'String', '3. Узлы Y (пусто = расчет f(x)):', ... % Метка
        'Position', [10 630 pW-25 20], 'HorizontalAlignment', 'left'); % Слева
    hEditY = uicontrol(panelControl, 'Style', 'edit', 'String', '', ... % Пусто по умолчанию
        'Position', [10 605 pW-25 25], 'BackgroundColor', 'white'); % Белый фон

    % 4. Секция ввода весов
    uicontrol(panelControl, 'Style', 'text', 'String', '4. Веса (4-я точка вес 25):', ... % Метка
        'Position', [10 575 pW-25 20], 'HorizontalAlignment', 'left', 'ForegroundColor', [0.6 0 0]); % Красный текст
    hEditW = uicontrol(panelControl, 'Style', 'edit', 'String', '1, 1, 1, 25, 1, 1, 1', ... % Веса
        'Position', [10 550 pW-25 25], 'BackgroundColor', 'white'); % Белый фон

    % --- ПАНЕЛЬ ИНТЕРПОЛЯЦИИ (Синяя группа) ---
    pInterp = uipanel(panelControl, 'Title', 'ИНТЕРПОЛЯЦИЯ', 'Units', 'pixels', ... % Подпанель
        'Position', [5 365 pW-15 175], 'BackgroundColor', [0.9 0.95 1]); % Синеватый фон
    hMenuI = uicontrol(pInterp, 'Style', 'popupmenu', ... % Список методов
        'String', {'Linear (Линейная)', 'Spline (Сплайн)', 'PCHIP (Эрмитова)', 'Lagrange (Лагранж)'}, ... % Методы
        'Position', [10 125 pW-40 25]); % Позиция
    uicontrol(pInterp, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ... % Кнопка
        'Position', [10 75 pW-40 40], 'BackgroundColor', [0.2 0.4 0.8], 'ForegroundColor', 'w', ... % Синий цвет
        'FontWeight', 'bold', 'Callback', @(~,~) mainCalculation('interp')); % Вызов расчета
    uicontrol(pInterp, 'Style', 'pushbutton', 'String', 'Таблица разностей', ... % Кнопка таблицы
        'Position', [10 25 pW-40 40], 'Callback', @(~,~) openDifferenceTable()); % Вызов таблицы

    % --- ПАНЕЛЬ АППРОКСИМАЦИИ (Зеленая группа) ---
    pApprox = uipanel(panelControl, 'Title', 'АППРОКСИМАЦИЯ', 'Units', 'pixels', ... % Подпанель
        'Position', [5 170 pW-15 190], 'BackgroundColor', [0.9 1 0.9]); % Зеленоватый фон
    hMenuA = uicontrol(pApprox, 'Style', 'popupmenu', ... % Список методов
        'String', {'МНК Полином', 'Взвешенный МНК', 'Экспонента', 'Степенная', 'Логарифм'}, ... % Список
        'Position', [10 145 pW-40 25]); % Позиция
    uicontrol(pApprox, 'Style', 'text', 'String', 'Степень m:', 'Position', [10 115 80 20], 'HorizontalAlignment', 'left'); % Текст m
    hEditM = uicontrol(pApprox, 'Style', 'edit', 'String', '2', 'Position', [95 110 50 25], 'BackgroundColor', 'w'); % Поле m
    uicontrol(pApprox, 'Style', 'pushbutton', 'String', 'РАССЧИТАТЬ', ... % Кнопка
        'Position', [10 60 pW-40 45], 'BackgroundColor', [0.1 0.6 0.1], 'ForegroundColor', 'w', ... % Зеленый цвет
        'FontWeight', 'bold', 'Callback', @(~,~) mainCalculation('approx')); % Вызов расчета
    uicontrol(pApprox, 'Style', 'pushbutton', 'String', 'Найти лучшую степень m', ... % Кнопка m
        'Position', [10 15 pW-40 35], 'Callback', @(~,~) findBestDegree()); % Вызов автоподбора

    % Поле для текстового вывода ошибок (RMSE)
    hResultText = uicontrol(panelControl, 'Style', 'edit', 'Max', 2, 'Min', 0, ... % Многострочное поле
        'Position', [10 10 pW-25 150], 'HorizontalAlignment', 'left', 'BackgroundColor', 'white'); % Снизу

    % Снова включаем адаптивность панели в процентах
    set(panelControl, 'Units', 'normalized'); % Чтобы растягивалось

    % Создаем область графика справа
    axesMain = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.33 0.1 0.64 0.85]); % Справа
    grid(axesMain, 'on'); % Координатная сетка

    % --- ФУНКЦИЯ СБОРА ДАННЫХ ИЗ ИНТЕРФЕЙСА ---
    function [x_data, y_data, w_data] = getInputsFromUI() % Начинаем сбор
        x_data = str2num(get(hEditX, 'String')); % Конвертируем строку X в массив
        y_input_str = get(hEditY, 'String'); % Читаем строку Y
        if isempty(y_input_str) % Если пользователь не ввел Y
            func_str = get(hEditFun, 'String'); % Берем строку f(x)
            f_handle = eval(['@(x)' func_str]); % Создаем функцию через eval
            y_data = f_handle(x_data); % Считаем значения функции в узлах
        else % Если пользователь ввел Y сам
            y_data = str2num(y_input_str); % Конвертируем строку Y в массив
        end % Конец проверки Y
        w_input_str = get(hEditW, 'String'); % Считываем строку весов
        if isempty(w_input_str) % Если весов нет
            w_data = ones(size(x_data)); % Все веса равны 1
        else % Если веса введены
            w_data = str2num(w_input_str); % Конвертируем в числа
        end % Конец проверки весов
    end % Конец функции

    % --- ГЛАВНАЯ ФУНКЦИЯ РАСЧЕТА ---
    function mainCalculation(mode) % Принимает 'interp' или 'approx'
        try % Начинаем защищенный блок
            [nodesX, nodesY, weights] = getInputsFromUI(); % Получаем все данные
            fineX = linspace(min(nodesX), max(nodesX), 500); % Создаем 500 точек для гладкости
            if strcmp(mode, 'interp') % Если ИНТЕРПОЛЯЦИЯ
                idxI = get(hMenuI, 'Value'); % Получаем номер метода из списка
                if idxI == 4 % Если выбран Лагранж
                    lagrangeCoeffs = polyfit(nodesX, nodesY, length(nodesX)-1); % Полином степени N-1
                    fineY = polyval(lagrangeCoeffs, fineX); % Считаем значения на сетке
                    methodName = 'Lagrange'; % Название для вывода
                else % Для других методов интерполяции
                    if idxI == 1, type = 'linear'; % 1-й индекс - Линейная
                    elseif idxI == 2, type = 'spline'; % 2-й индекс - Сплайн
                    else, type = 'pchip'; end % 3-й индекс - Эрмитова
                    fineY = interp1(nodesX, nodesY, fineX, type); % Встроенная функция interp1
                    methodName = type; % Название для вывода
                end % Конец выбора метода
            else % Если АППРОКСИМАЦИЯ
                m = str2double(get(hEditM, 'String')); % Считываем m
                idxA = get(hMenuA, 'Value'); % Получаем номер метода
                % ФИКС НУЛЕЙ (epsilon-shift)
                s = 0.05; % Маленькое смещение
                safeX = nodesX; % Копия X
                safeY = nodesY; % Копия Y
                safeFineX = fineX; % Копия сетки
                if idxA >= 3 % Если это Экспонента, Степень или Лог
                    if any(nodesX <= 0), safeX = nodesX + s; safeFineX = fineX + s; end % Сдвиг X
                    if any(nodesY <= 0), safeY = nodesY + s; end % Сдвиг Y
                end % Конец фикса
                switch idxA % Выбираем математику
                    case 1 % МНК Полином
                        p = polyfit(nodesX, nodesY, m); % Встроенный МНК
                        fineY = polyval(p, fineX); % Итоговая линия
                        methodName = 'МНК Полином'; % Имя
                    case 2 % Взвешенный МНК
                        V = zeros(length(nodesX), m+1); % Матрица Вандермонда
                        for i = 0:m % Цикл по степеням
                            V(:, i+1) = nodesX(:).^i; % Наполняем x^0, x^1...
                        end % Конец матрицы
                        coeffsW = lscov(V, nodesY', weights); % Решение с весами (WLS)
                        fineY = (zeros(500, m+1) + fineX(:).^(0:m)) * coeffsW; % Линия
                        methodName = 'Взвешенный МНК'; % Имя
                    case 3 % Экспонента
                        p_exp = polyfit(safeX, log(abs(safeY)), 1); % Линеаризация ln(y)
                        fineY = exp(p_exp(2)) * exp(p_exp(1) * safeFineX); % Перевод обратно
                        methodName = 'Экспонента'; % Имя
                    case 4 % Степенная
                        p_pow = polyfit(log(abs(safeX)), log(abs(safeY)), 1); % ln(y) по ln(x)
                        fineY = exp(p_pow(2)) * safeFineX.^p_pow(1); % Перевод обратно
                        methodName = 'Степенная'; % Имя
                    case 5 % Логарифм
                        p_log = polyfit(log(abs(safeX)), safeY, 1); % y по ln(x)
                        fineY = p_log(1) * log(abs(safeFineX)) + p_log(2); % Линия логарифма
                        methodName = 'Логарифм'; % Имя
                end % Конец switch
            end % Конец режима
            cla(axesMain); % Очищаем график перед новой отрисовкой
            plot(axesMain, fineX, fineY, 'r', 'LineWidth', 2); % Рисуем красную линию
            hold(axesMain, 'on'); % Режим наложения
            scatter(axesMain, nodesX, nodesY, 80, 'b', 'filled'); % Рисуем синие точки узлов
            grid on; title(axesMain, ['Результат: ' methodName]); % Сетка и титул
            % Считаем погрешность в узлах
            yCalcNodes = interp1(fineX, fineY, nodesX, 'linear', 'extrap'); % Снимаем значения
            rmse = sqrt(mean((nodesY - yCalcNodes).^2)); % Формула RMSE
            maxE = max(abs(nodesY - yCalcNodes)); % Формула Max Error
            outputStr = sprintf('МЕТОД: %s\nRMSE: %.4e\nMAX ERR: %.4e', methodName, rmse, maxE); % Строка текста
            set(hResultText, 'String', outputStr); % Выводим текст в поле UI
        catch ME % Если случилась ошибка
            errordlg(['Ошибка в расчете: ' ME.message]); % Показываем окно ошибки
        end % Конец блока try
    end % Конец главной функции

    function openDifferenceTable() % Функция таблицы разностей
        [~, nY, ~] = getInputsFromUI(); % Получаем значения Y
        nP = length(nY); % Количество точек
        D = zeros(nP, nP); % Пустая матрица
        D(:, 1) = nY(:); % 1-й столбец это Y
        for j = 2:nP % Цикл по порядку разности
            for i = 1:nP - j + 1 % Цикл по строкам
                D(i, j) = D(i+1, j-1) - D(i, j-1); % Разность соседних элементов
            end % Конец строк
        end % Конец столбцов
        fT = figure('Name', 'Таблица разностей', 'Position', [400 400 600 450]); % Окно
        uitable(fT, 'Data', D, 'Units', 'normalized', 'Position', [0 0 1 1]); % Отрисовка
    end % Конец функции разностей

    function findBestDegree() % Функция подбора m
        [nX, nY, ~] = getInputsFromUI(); % Узлы X и Y
        bestM = 1; % По умолчанию
        minR = inf; % Ошибка бесконечность
        for m_test = 1:length(nX)-1 % Перебор степеней
            coeffs_test = polyfit(nX, nY, m_test); % Считаем МНК
            y_test = polyval(coeffs_test, nX); % Значения в узлах
            r_test = sqrt(mean((nY - y_test).^2)); % RMSE
            if r_test < minR % Если ошибка меньше
                minR = r_test; bestM = m_test; % Запоминаем
            end % Конец условия
        end % Конец цикла
        set(hEditM, 'String', num2str(bestM)); % Обновляем интерфейс
    end % Конец функции подбора m
end % Конец всего файла App_BuiltIn