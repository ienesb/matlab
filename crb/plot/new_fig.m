function fig = new_fig(pos)
% NEW_FIG  Invisible white figure with the light theme forced (batch-mode
% MATLAB R2025 otherwise renders dark axes on a white figure background).
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', pos);
    try, fig.Theme = 'light'; catch, end                       %#ok<TRYNC>  (R2025a+)
end
