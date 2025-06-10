%% Settings for publication quality figures
% interpret text with LaTeX. Can find others with `get(groot,'factory')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); % most plotting uses this
set(groot, 'defaultTextInterpreter', 'latex'); 
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultGraphplotInterpreter','latex');
set(groot, 'defaultConstantlineInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextboxshapeInterpreter','latex');
% There are more interpreters I must set to latex (tiled plot?)

set(groot, 'defaultTiledlayoutTileSpacing', 'tight');
set(groot, 'defaultTiledlayoutPadding', 'tight');

% line settings
colororder({'red','green','blue','cyan','magenta'}) % don't use yellow for the line color
set(groot, 'defaultLineLineWidth', 2)

% default title, ticks, xlabel, ylabel, and legend size
set(groot, 'defaultAxesFontSize',14);
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultConstantlineLineWidth', 2);
close all

% symbolic toolbox preferences
% sympref('PolynomialDisplayStyle','descend');
sympref('AbbreviateOutput',false);


% these legend setting are seemingly useless
% set(groot, 'defaultLegendFontSize', 14);
% set(groot, 'defaultLegendInterpreter', 'latex');
% Also, we can't set font weight when using the latex interpreter
% set(groot,'defaultTextFontWeight','bold')