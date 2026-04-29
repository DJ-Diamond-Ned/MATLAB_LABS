clear; clc; close all;

% 1. Входные данные (6 точек)
T_points = [45, 165, 315, 465, 615, 830]; % °F
T_full = linspace(45, 830, 24); % сетка для 24 гипокомпонентов

% Данные из таблиц (6 значений в каждой строке)
groups(1).name = 'NaphthenesHydrogenation';
groups(1).compNames = {'Cyclopentanes', 'Cyclohexanes', 'NaphtheneAromatic', 'DiNaphtheneAromatic', 'TriNaphtheneAromatic'};
groups(1).data = [
    1.083E-10, 1.083E-10, 9.915E-01, 5.016E-01, 2.602E-01, 1.808E-01;
    9.972E-01, 9.972E-01, 8.458E-03, 4.926E-01, 3.045E-01, 2.227E-01;
    2.433E-03, 2.433E-03, 4.268E-12, 5.723E-03, 3.597E-01, 4.299E-01;
    3.295E-04, 3.295E-04, 4.268E-12, 6.489E-05, 7.195E-02, 1.444E-01;
    2.588E-06, 2.588E-06, 4.268E-12, 5.096E-07, 3.628E-03, 2.218E-02];

groups(2).name = 'AromHydrog';
groups(2).compNames = {'MonoAromatic', 'DiAromatic', 'TriAromatic', 'TetraAromatic'};
groups(2).data = [
    9.916E-01, 9.916E-01, 1.000E+00, 9.924E-01, 5.367E-01, 3.057E-01;
    8.368E-03, 8.368E-03, 2.069E-10, 7.621E-03, 4.181E-01, 3.892E-01;
    3.724E-10, 3.724E-10, 8.311E-06, 2.168E-13, 4.523E-02, 2.340E-01;
    3.724E-10, 3.724E-10, 1.589E-08, 2.168E-13, 1.732E-13, 7.107E-02];

groups(3).name = 'HDN';
groups(3).compNames = {'Pyrroles', 'Pyridines', 'Indoles', 'Quinolones', 'Carbazoles'};
groups(3).data = [
    4.480E-05, 2.763E-01, 4.480E-05, 6.729E-01, 5.863E-01, 3.690E-01;
    4.480E-05, 7.226E-01, 4.480E-05, 3.271E-01, 2.481E-01, 1.583E-01;
    9.998E-01, 5.414E-08, 9.998E-01, 4.374E-12, 1.276E-01, 3.584E-01;
    4.480E-05, 1.057E-03, 4.480E-05, 1.118E-05, 3.801E-02, 1.009E-01;
    7.867E-05, 5.414E-08, 7.867E-05, 4.374E-12, 2.259E-06, 1.349E-02];

groups(4).name = 'HDS';
groups(4).compNames = {'Sulfides', 'Thiophenes', 'Benzothiophene', 'Dibenzothiophene'};
groups(4).data = [
    2.062E-06, 5.742E-01, 2.062E-06, 3.203E-01, 2.065E-01, 1.045E-01;
    2.062E-06, 4.258E-01, 2.062E-06, 6.788E-01, 4.099E-01, 2.102E-01;
    9.982E-01, 5.820E-08, 9.982E-01, 8.592E-04, 3.786E-01, 4.722E-01;
    1.823E-03, 5.820E-08, 1.823E-03, 3.657E-12, 4.941E-03, 2.131E-01];

% 2. Расчет и нормировка
for g = 1:length(groups)
    numComp = length(groups(g).compNames);
    interp_vals = zeros(numComp, length(T_full));
    
    figure('Name', ['Группа: ' groups(g).name], 'Color', 'w');
    hold on; grid on; box on;
    colors = lines(numComp);
    
    for c = 1:numComp
        y = groups(g).data(c, :);
        y_raw = interp1(T_points, y, T_full, 'pchip');
        y_raw(y_raw < 0) = 0;
        interp_vals(c, :) = y_raw;
    end
    
    col_sums = sum(interp_vals, 1);
    final_vals = bsxfun(@rdivide, interp_vals, col_sums);
    
    % Отрисовка
    for c = 1:numComp
        plot(T_full, final_vals(c, :), 'Color', colors(c,:), 'LineWidth', 2, 'DisplayName', groups(g).compNames{c});
        scatter(T_points, groups(g).data(c,:), 50, colors(c,:), 'filled', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    title(['Состав группы: ', groups(g).name]);
    xlabel('T, ^oF'); ylabel('Мольная доля');
    legend('Location', 'bestoutside');
    
    % Вывод таблицы в командное окно
    fprintf('\n--- Группа: %s ---\n', groups(g).name);
    VarNames = cell(1, 24);
    for i=1:24, VarNames{i} = sprintf('T_%.0f', T_full(i)); end
    Tbl = array2table(final_vals, 'RowNames', groups(g).compNames, 'VariableNames', VarNames);
    disp(Tbl);
end