% combined_corr_curves.m
% In this script, six values are closely approximated via golden section search
% - α value for which correlation coefficient ρ is strongest between E_π and T_{a}^1
% - α value for which correlation coefficient ρ is strongest between E_π and T_{a}^2
% E_π = Pi Electronic Energy
% T_{a}^1 = First General Temperature Indices
% T_{a}^2 = Second General Temperature Indices
% Additionally, curves ρ against α near these values are plotted in 2 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = true; % Set to true to auto-save plots

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

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
% Usage: getIndexFns{n}(alpha, nu_exp) | n=1:T_{a}^1 , n=2:T_{a}^2, a = alpha
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (nu_exp - 2)), ((5 .* nu_exp - 12) ./ ((nu_exp - 2) .* (nu_exp - 3))), (6 ./ (nu_exp - 3))].^alpha, 2); % First General Temperature Indices (T_{a}^1)
    @(alpha) sum(d_f .* [(4 ./ ((nu_exp - 2).^2)), (6 ./ ((nu_exp - 2) .* (nu_exp - 3))), (9 ./ ((nu_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_{a}^2)
}';

% Cell containing their respective labels
indexName = {"T^1", "T^2"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % one (E_{π})
numIndices = size(getIndexFns,2); % two (T_{a}^1  & T_{a}^2)
numCurves = numData*numIndices;   % two (2 curves)

% All x in visible ranges (both plots - near and far)
xnear = linspace(-0.23, 0.2, 900);
xfar = linspace(-20,20,800); % for E_{π}

% Do the same procedure for each experimental data i.e., bp
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % WARNING: these xmeet1-ymeet2 values are hardcoded, computed separately
  %          ... These are coordinates where ρ-α of E_{π}-T_{a}^1
  %          ... intersects ρ-α of E_{π}-T_{a}^2
  xmeet1 = -0.05966;  % for E_{π}
  ymeet1 = 0.99715;   % for E_{π}
  xmeet2 = -0.0005;  % for E_{π}
  ymeet2 = 0.99691;   % for E_{π}

  ybox = [
  0.99750; % for figure 2 (blue dotted line)
  ](ii);
  % Plot the blue dashed box (before drawing the curves so it appear beneath)
  plot([xmeet1 xmeet1 0 0], [0 ybox ybox 0], '--b', 'LineWidth', lineWidth);


  yend = 0; % <-- to be assigned some value later for adjusting visible range

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between E_π with specified α
    %                                and T_{a}^1 /T_{a}^2 (depending on n)
    get_indices_vals = @(alpha) getIndexFns{n}(alpha)(!isnan(expData{1,ii}));
    ccFn = @(alpha) corrcoef(
      get_indices_vals(alpha), % Either T_{a}^1  or T_{a}^2
      expData{1,ii}(!isnan(expData{1,ii})) % E_{π}
    )(1,2);

    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear);
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(GoldenSectionSearch_Maximum(ccFn, xnear(1), xnear(end), 1e-15));
    peakCorrCoeff = ccFn(peakAlpha);

    % Display peak alpha and peak correlation coefficient in the command window
    disp(['Peak Alpha: ', num2str(peakAlpha)]);
    disp(['Peak Correlation Coefficient: ', num2str(peakCorrCoeff)]);


    % Generate curve label [E_{π}] [T_{a}^1/T_{a}^2]
    curveLabels{n} = sprintf("%s and %s_a", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear, ynear, '-', 'LineWidth', lineWidth);

    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)

     yend = max(yend, ynear(end)); % y value to be used as visible y lower bound
  end

  % Mark and write on the plot the limits of alpha where T_alpha_1  is better than T_α^2
  plot([xmeet1 xmeet2], [ymeet1 ymeet2], '*b',
       'MarkerSize', 16, 'LineWidth', lineWidth/1.5); % Mark with blue asterisks
  text(xmeet1, ymeet1, {'', sprintf(" (−%s, %s)", as_4_dp_str(abs(xmeet1)), as_4_dp_str(ymeet1))},
       'VerticalAlignment', 'top', 'Color', [0, 0, 0.8]); % Write blue text
  text(xmeet2, ymeet2, {'', sprintf("(0, %s) ", as_4_dp_str(ymeet2))},
       'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'Color', [0, 0, 0.8]);

  % Label this expData's zoomed-in plot
  xlabel('α');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "northeast"); % the location of the legend box

  ybox_space = 0.000005; % spacing between upper y-value with the blue dotted lines (figure 3)
  axis([xnear(1) xnear(end) yend ybox+ybox_space]); % Enforce figure's visible range
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('α');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  set(leg, 'location', "southeast");

  hold off;
end

for ii = 1:2
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));

  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_E_{π}_indices_FAR.png");
  saveas(figure(2), "01_comb_ccurves_E_{π}_indices_NEAR.png");
end
