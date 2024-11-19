% scatter_plots.m
% In this script, scatter plots are generated for the pair
%  - E_{π} and T_{a}^1_{-0.1589}
%  - E_{π} and T_{a}^2_{-0.0778}
% ... with their respective regression lines.

close all; % Before drawing, close any figures already opened
clear;     % Clear all variables

% CONFIG: Line width and font size for each curve in drawn figures
lineWidth = 2;
fontSize = 26;
% Save plots to images? Set to true if yes
saveToFile = false;

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.06f', round(x*(10^6))/(10^6));

% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([% Pi Electroni Energy
   8.0000  13.6832 19.3137 19.4483 24.9308 25.1875 25.1012 25.1922 25.2745 22.5055
   30.5440 30.7255 30.8805 30.8795 30.7627 30.9990 30.9362 30.9386 30.9432 30.8390
   30.9418 30.8336 28.2453 28.3361 28.2220 36.6814 31.4251 36.1557 34.5718 46.4974
   ]', 30, 1),
   "E_{π}" % Their labels
 };

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
 d_f = [
   6 6 6 7 6  8 7  8 9 6 6  7  8  8  7  10 9  8  9  8  9  8  8 8 7  6  7  6  6  4
   0 4 8 6 12 8 10 8 6 8 16 14 12 12 14 8  10 12 10 12 10 10 8 8 10 12 10 20 12 16
   0 1 2 3 3  5 4  5 6 5 4  5  6  6  5  8  7  6  7  6  7  6  8 8 7  12 10 5  12 19
 ]'; % Used for computing indices based on edge endpoint degree partitions

% Define nu_exp here before using it in getIndexFns Lower 30 Bhs (nu_exp: number of vertices)
nu_exp = [
    6 10 14 14 18 18 18 18 18 18 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 24 22 26 24 32
]';

% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(alpha, nu_exp) | n=1:T_{a}^1, n=2:T_{a}^2, a = alpha
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (nu_exp - 2)), ((5 .* nu_exp - 12) ./ ((nu_exp - 2) .* (nu_exp - 3))), (6 ./ (nu_exp - 3))].^alpha, 2); % First General Temperature Indices (T_{a}^1)
    @(alpha) sum(d_f .* [(4 ./ ((nu_exp - 2).^2)), (6 ./ ((nu_exp - 2) .* (nu_exp - 3))), (9 ./ ((nu_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_{a}^2)
}';

% Cell containing their respective labels
indexName = {"T_1" "T_2"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % one (E_{π})
numIndices = size(getIndexFns,2); % two (T_{a}^1 & T_{a}^2)
numCurves = numData*numIndices;   % two (2 curves)

for edn = 1:numData % edn = experimental data number | 1=E_{π}
  for fnn = 1:numIndices % fnn = function number | n=1:T_{a}^1, n=2:T_{a}^2
    ccFn = @(alpha) corrcoef( % Gets corrcoef between bp and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);         % For the regression line
    x = [getIndexFns{fnn}(peakAlpha)(1) max(getIndexFns{fnn}(peakAlpha))];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth', lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s^{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  % Save each figure to a separate file
  saveas(figure(1), "Scatter_E_T^1_a.png");
  saveas(figure(2), "Scatter_E_T^2_a.png");
end
