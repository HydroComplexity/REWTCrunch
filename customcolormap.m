function [customCMap] = customcolormap(data,white,datlargest,datsmallest)
indexValue = white;             % value to be represented with WHITE
topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])

% Calculate where proportionally indexValue lies between minimum and
% maximum values
L = length(data(:));
largest = datlargest; %max(max(data));
smallest = datsmallest; %min(min(data));
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100.*index)',...
    linspace(bottomcolor(2),indexColor(2),100.*index)',...
    linspace(bottomcolor(3),indexColor(3),100.*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
    linspace(indexColor(2),topColor(2),100*(L-index))',...
    linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps
end