clear; clc; close all;

% 1. ИСХОДНЫЕ ДАННЫЕ
ObservedTemperatures = [45, 165, 315, 465, 615, 830]; 

Temperatures_24_Components = linspace(45, 830, 24); 

% Сетка из 500 точек для построения красивого гладкого графика
Temperatures_Smooth_Graph = linspace(45, 830, 500); 

% Структура с данными по каждой реакционной группе
ReactionGroups(1).name = 'NaphthenesHydrogenation';
ReactionGroups(1).compNames = {'Cyclopentanes', 'Cyclohexanes', 'NaphtheneAromatic', 'DiNaphtheneAromatic', 'TriNaphtheneAromatic'};
ReactionGroups(1).data = [1.083E-10, 1.083E-10, 9.915E-01, 5.016E-01, 2.602E-01, 1.808E-01; 9.972E-01, 9.972E-01, 8.458E-03, 4.926E-01, 3.045E-01, 2.227E-01; 2.433E-03, 2.433E-03, 4.268E-12, 5.723E-03, 3.597E-01, 4.299E-01; 3.295E-04, 3.295E-04, 4.268E-12, 6.489E-05, 7.195E-02, 1.444E-01; 2.588E-06, 2.588E-06, 4.268E-12, 5.096E-07, 3.628E-03, 2.218E-02];

ReactionGroups(2).name = 'AromHydrog';
ReactionGroups(2).compNames = {'MonoAromatic', 'DiAromatic', 'TriAromatic', 'TetraAromatic'};
ReactionGroups(2).data = [9.916E-01, 9.916E-01, 1.000E+00, 9.924E-01, 5.367E-01, 3.057E-01; 8.368E-03, 8.368E-03, 2.069E-10, 7.621E-03, 4.181E-01, 3.892E-01; 3.724E-10, 3.724E-10, 8.311E-06, 2.168E-13, 4.523E-02, 2.340E-01; 3.724E-10, 3.724E-10, 1.589E-08, 2.168E-13, 1.732E-13, 7.107E-02];

ReactionGroups(3).name = 'HDN';
ReactionGroups(3).compNames = {'Pyrroles', 'Pyridines', 'Indoles', 'Quinolones', 'Carbazoles'};
ReactionGroups(3).data = [4.480E-05, 2.763E-01, 4.480E-05, 6.729E-01, 5.863E-01, 3.690E-01; 4.480E-05, 7.226E-01, 4.480E-05, 3.271E-01, 2.481E-01, 1.583E-01; 9.998E-01, 5.414E-08, 9.998E-01, 4.374E-12, 1.276E-01, 3.584E-01; 4.480E-05, 1.057E-03, 4.480E-05, 1.118E-05, 3.801E-02, 1.009E-01; 7.867E-05, 5.414E-08, 7.867E-05, 4.374E-12, 2.259E-06, 1.349E-02];

ReactionGroups(4).name = 'HDS';
ReactionGroups(4).compNames = {'Sulfides', 'Thiophenes', 'Benzothiophene', 'Dibenzothiophene'};
ReactionGroups(4).data = [2.062E-06, 5.742E-01, 2.062E-06, 3.203E-01, 2.065E-01, 1.045E-01; 2.062E-06, 4.258E-01, 2.062E-06, 6.788E-01, 4.099E-01, 2.102E-01; 9.982E-01, 5.820E-08, 9.982E-01, 8.592E-04, 3.786E-01, 4.722E-01; 1.823E-03, 5.820E-08, 1.823E-03, 3.657E-12, 4.941E-03, 2.131E-01];

% 2. ОСНОВНОЙ РАСЧЕТНЫЙ ЦИКЛ
for g = 1:length(ReactionGroups)
    CurrentGroupCount = length(ReactionGroups(g).compNames);
    
    % Создаем матрицы для хранения результатов интерполяции до нормировки
    RawFractions_Table_24 = zeros(CurrentGroupCount, length(Temperatures_24_Components));
    RawFractions_Smooth_Graph = zeros(CurrentGroupCount, length(Temperatures_Smooth_Graph));
    
    for c = 1:CurrentGroupCount
        % Берем данные Unisim для конкретного вещества
        ObservedMolarFractions = ReactionGroups(g).data(c, :);
        
        % ПРИМЕНЯЕМ ИНТЕРПОЛЯЦИЮ PCHIP (степень <= 3)
        CalculatedFractions_24 = interp1(ObservedTemperatures, ObservedMolarFractions, Temperatures_24_Components, 'pchip');
        CalculatedFractions_24(CalculatedFractions_24 < 0) = 0; % Убираем отрицательные значения
        RawFractions_Table_24(c, :) = CalculatedFractions_24;
        
        % 2. Расчет для плавного графика (500 точек)
        CalculatedFractions_Smooth = interp1(ObservedTemperatures, ObservedMolarFractions, Temperatures_Smooth_Graph, 'pchip');
        CalculatedFractions_Smooth(CalculatedFractions_Smooth < 0) = 0;
        RawFractions_Smooth_Graph(c, :) = CalculatedFractions_Smooth;
    end
    
    % Нормировка для таблицы
    SumOfAllComponents_Table = sum(RawFractions_Table_24, 1);
    Fractions_Normalized_Table = bsxfun(@rdivide, RawFractions_Table_24, SumOfAllComponents_Table);
    
    % Нормировка для графика
    SumOfAllComponents_Graph = sum(RawFractions_Smooth_Graph, 1);
    Fractions_Normalized_Smooth = bsxfun(@rdivide, RawFractions_Smooth_Graph, SumOfAllComponents_Graph);
    
    % 3. ВИЗУАЛИЗАЦИЯ И ВЫВОД
    
    % Построение графика
    figure('Name', ['График группы: ' ReactionGroups(g).name], 'Color', 'w');
    hold on; grid on; box on;
    ColorPalette = lines(CurrentGroupCount);
    
    for c = 1:CurrentGroupCount
        % Отрисовка плавной линии по 500 точкам
        plot(Temperatures_Smooth_Graph, Fractions_Normalized_Smooth(c, :), 'Color', ColorPalette(c,:), 'LineWidth', 2, 'DisplayName', ReactionGroups(g).compNames{c});
        % Нанесение исходных 6 точек (маркеров)
        scatter(ObservedTemperatures, ReactionGroups(g).data(c,:), 60, ColorPalette(c,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    title(['Фракционный состав группы: ', ReactionGroups(g).name]);
    xlabel('Температура кипения T, ^oF'); ylabel('Мольная доля');
    legend('Location', 'bestoutside');
    xlim([0 900]); ylim([0 1.1]);
    
    % Формирование таблицы в консоли
    fprintf('\n--- ИТОГОВАЯ ТАБЛИЦА: ГРУППА %s ---\n', ReactionGroups(g).name);
    ColumnTitles = cell(1, 24);
    for i=1:24, ColumnTitles{i} = sprintf('T_%.0f', Temperatures_24_Components(i)); end
    
    FinalTable = array2table(Fractions_Normalized_Table, 'RowNames', ReactionGroups(g).compNames, 'VariableNames', ColumnTitles);
    disp(FinalTable);
end