% good_alpha_intervals.m
% In this script, 2 plots are generated which shows α interals for
% - good correlation coefficient ρ between E_{π} and T_{a}^1
% - good correlation coefficient ρ between E_{π} and T_{a}^2
% The limits of the intervals are also printed to console

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 3;
fontSize = 24;
saveToFile = false; % Set to true to auto - save plots
% Note: Dimensions of resulting images are scaled according to each window size.
%       To control the dimensions, after running the whole script, resize each
%       ... figure window and then run only the saveas functions
%       ... manually, by selection, at the end of this script

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 22 lower benzenoids
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
% Usage: getIndexFns{n}(alpha, nu_exp) | n=1:T_{a}^1, n=2:T_{a}^2
getIndexFns = {
    @(alpha) sum(d_f .* [(4 ./ (nu_exp - 2)), ((5 .* nu_exp - 12) ./ ((nu_exp - 2) .* (nu_exp - 3))), (6 ./ (nu_exp - 3))].^alpha, 2); % First General Temperature Indices (T_{a}^1)
    @(alpha) sum(d_f .* [(4 ./ ((nu_exp - 2).^2)), (6 ./ ((nu_exp - 2) .* (nu_exp - 3))), (9 ./ ((nu_exp - 3).^2))].^alpha, 2); % Second General Temperature Indices (T_{a}^2)
}';

% Cell containing their respective labels
indexName = {"T^t_1", "T^t_2"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % one (E_{π})
numIndices = size(getIndexFns,2); % two (T_{a}^1 & T_{a}^2)
numCurves = numData*numIndices;   % two (2 curves)

% Boundaries for visible intervals, for each index-property pair
%             T_{a}^1  T_{a}^2
xstart = [   -2     -0.8];   % E_{π}
xend =   [    0.8     0.4];   % E_{π}
ystart = [   1   0.998]  % E_{π}
yend =   [   0.94    0.96]; % E_{π}

% Exact rho value considered good for
%            E_{π}
a_goodrho = [0.98];

% Colors (different shades of cyan and green)
colShaded = {[0.85, 1, 1]; [0.85, 1, 0.85]};
colIndicatorVert = {[0.2, 0.55, 0.55]; [0, 0.5, 0]};
colIndicatorHorz = {[0.35, 0.75, 0.75]; [0.45, 0.7, 0.45]};
colCurve = {[0, 0.75, 0.75]; [0, 0.5, 0]};

  % Do the same procedure for each experimental data i.e., E_{π}
for ii = 1:numData
  % Create a new figure for each experimental data
  figure('Name', sprintf('Correlation between %s and %s', indexName{1}, expData{1,ii})); % Adjusted line
  hold on;

   % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between E_{π} with specified α
    %                          and T_{a}^1/T_{a}^2 (depending on n)
    ccFn = @(alpha) corrcoef(
      getIndexFns{n}(alpha)(~isnan(expData{1,ii})),
      expData{1,ii}(~isnan(expData{1,ii}))
    )(1,2);

    this_fig = figure((ii-1)*numIndices+n);
    hold on;

    % Get Interval limits and print them to console

    % Get alpha for highest rho first
    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -4, 0, 1e-15));

    % Prepare a funcion to calc |rho(a)-goodrho|
    ccFn_good = @(a)(-abs(ccFn(a)-a_goodrho(ii)));

    % and func to get the limit, i.e., value of alpha where rho is 0.98 or 0.99
    getLimitFromInterval = @(lb, ub) mean(
      GoldenSectionSearch_Maximum(ccFn_good, lb, ub, 1e-15));

    a_lb = getLimitFromInterval(peakAlpha-1, peakAlpha); % Search to the left
    a_ub = getLimitFromInterval(peakAlpha, peakAlpha+1); % Search to the right

    % Write the intervals in console
    disp(sprintf("ρ(%s,%s_α) ≥ %.02f when α ∈ [%.08f, %.08f]",
         expData{2,ii}, indexName{n}, a_goodrho(ii), a_lb, a_ub));

    % Plot the actual curve, but exclude good alpha range (<- separately drawn)
    % generate x values
    x = [linspace(xstart(ii,n),a_lb,300), linspace(a_ub,xend(ii,n),300)];
    % generate the corresponding y values
    y = arrayfun(ccFn, x);
    plot(x, y, '-', 'LineWidth', lineWidth);

    % Shade the area inside the good alpha interval
    u0 = a_lb;         u_width = a_ub-a_lb;
    v0 = ystart(ii,n); v_height = yend(ii,n) - ystart(ii,n);
    rectangle('Position', [u0, v0, u_width, v_height], 'FaceColor', colShaded{ii}, 'LineStyle', 'none');

    % Draw the indicator lines for the good alpha interval
    % Vertical dashed lines:
    plot([a_lb a_lb], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
    plot([a_ub a_ub], [ystart(ii,n) yend(ii,n)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});

    % Horizontal black dashed line, horizontal colored dashed line:
    plot([xstart(ii,n) a_lb], [a_goodrho(ii) a_goodrho(ii)],
         '--k', 'LineWidth', lineWidth/1.75);
    plot([a_lb a_ub], [a_goodrho(ii) a_goodrho(ii)],
         '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});

    % Write on the plot the interval limits
    text(a_lb, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_lb)))}, 'VerticalAlignment', 'top', 'Rotation', 90);
    text(a_ub, yend(ii,n), {'', sprintf("  α=−%s", as_4_dp_str(abs(a_ub)))}, 'VerticalAlignment', 'top', 'Rotation', 90);

    % Finally, plot the colored curve within the good alpha range
    x_in = linspace(a_lb,a_ub,400);
    y_in = arrayfun(ccFn, x_in);
    plot(x_in, y_in, '-', 'LineWidth', lineWidth, 'Color', colCurve{ii});

    % Also highlight the good interval on the x-axis
    plot([a_lb, a_ub], [yend(ii,n), yend(ii,n)], '-',
         'LineWidth', lineWidth, 'Color', colCurve{ii});

    % Label the plot                       [E_{π} ]        [T_{a}^1,T_{a}^2]
    title(sprintf('between %s and %s', expData{2,ii}, indexName{n}));
    xlabel('θ');
    ylabel('ρ');
    drawnow;
    axis([xstart(ii,n) xend(ii,n) yend(ii,n) ystart(ii,n)]);


    % Replace all hypens with minuses
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    % Set the font size
    set(findall(this_fig,'-property','FontSize'),'FontSize', fontSize);
    drawnow;

    hold off;
  end
end

if saveToFile
  saveas(figure(1), "02_good_a_intervals_E_T^1_a.png");
  saveas(figure(2), "02_good_a_intervals_E_T^2_a.png");
  saveas(figure(3), "02_good_a_intervals_E_SO_a.png");
end
