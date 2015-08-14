function subplots(data, layout, names, formatFn)

% Plot each element of a cell array or struct in its own subplot.
%
% subplots(data, layout, names, formatFn)
%
% formatFn will be called for each subplot and it can do whatever it wants
% (e.g. set axes, label axes) with three arguments like subplot:
%   formatFn(row, col, i)

if ~exist('layout', 'var'), layout = []; end
if ~exist('names', 'var'), names = {}; end
if ~exist('formatFn', 'var') || isempty(formatFn), formatFn = @(r,c,i) []; end

if isNoDisplay()
    % Don't bother...
    return
end

if isstruct(data)
    names = fieldnames(data);
    data  = struct2cell(data);
elseif ~iscell(data)
    data = {data};
end

nPlots = length(data);
if length(layout) < 2 || all(layout < 0)
    nRows = round(sqrt(nPlots));
    nCols = ceil(nPlots / nRows);
elseif layout(1) < 0
    nCols = layout(2);
    nRows = ceil(nPlots / nCols);
elseif layout(2) < 0
    nRows = layout(1);
    nCols = ceil(nPlots / nRows);
else
    nRows = layout(1);
    nCols = layout(2);
end
if isempty(names)
    names = cellfun(@(x) num2str(x), num2cell(1:nPlots), 'UniformOutput', false);
end

clf
for i = 1:nPlots
    subplot(nRows, nCols, i);
    x = data{i};
    if size(x, 2) < 5
        plot(x)
        axis tight
    elseif size(x, 1) < 5
        plot(x')
        axis tight
    else
        imagesc(x)
        axis xy
    end
    colorbar
    title(names{i}, 'Interpreter', 'none');
    row = floor((i-1) / nCols) + 1;
    col = mod(i-1, nCols) + 1;
    formatFn(row, col, i);
end

try
    %subplot 111
catch
    warning('Could not restore subplot to 111')
end
