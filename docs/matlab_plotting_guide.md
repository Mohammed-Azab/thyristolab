# MATLAB Plotting Guide for Research Papers

This guide documents the plotting conventions and syntax used in the ThyristoLab project for generating publication-quality figures.

## Table of Contents
- [Figure Setup](#figure-setup)
- [LaTeX Formatting](#latex-formatting)
- [Line Styles and Markers](#line-styles-and-markers)
- [Axis Configuration](#axis-configuration)
- [Legends](#legends)
- [Color Schemes](#color-schemes)
- [Scale Types](#scale-types)
- [Best Practices](#best-practices)

---

## Figure Setup

### Basic Figure Creation
```matlab
% Create a figure with specified size
figure('Name', 'Figure Title', 'Position', [100, 100, 800, 600]);
% Position: [left, bottom, width, height] in pixels
```

### Common Figure Sizes
- **Standard plot**: `[100, 100, 800, 600]`
- **Wide plot**: `[100, 100, 900, 600]`
- **Large multi-panel**: `[100, 100, 1200, 800]`

### Tiled Layouts (Multiple Subplots)
```matlab
% Create a 2x2 tiled layout
figure('Position', [100, 100, 1200, 800]);
t_layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_layout, 'Overall Title', 'Interpreter', 'latex', 'FontSize', 18);

% Add plots to tiles
nexttile;
plot(x, y, 'b-', 'LineWidth', 2);
% ... configure plot ...

nexttile;
% ... next plot ...
```

---

## LaTeX Formatting

### Enabling LaTeX Interpreter
Always use LaTeX for mathematical notation in research papers:

```matlab
% Axis labels
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Power Loss (W)', 'Interpreter', 'latex', 'FontSize', 16);

% Title
title('Power Losses vs Firing Angle', 'Interpreter', 'latex', 'FontSize', 18);

% Legend
legend('$P_{\mathrm{batt}}$', '$P_{\mathrm{total}}$', 'Interpreter', 'latex', 'FontSize', 14);
```

### Common LaTeX Symbols
- Greek letters: `$\alpha$`, `$\beta$`, `$\theta$`, `$\omega$`
- Subscripts: `$V_{\mathrm{out}}$`, `$I_{\mathrm{rms}}$`, `$R_{\mathrm{bat}}$`
- Superscripts: `$I^2$`, `$e^{-t}$`, `$x^{\circ}$` (degree symbol)
- Fractions: `$\frac{V_m}{2\pi}$`
- Mathematical operators: `$\times$`, `$\div$`, `$\pm$`

### Formatting Text with Math
Use `\mathrm{}` for non-italic text in equations:
```matlab
'$P_{\mathrm{batt}} = I_{\mathrm{rms}}^2 R_{\mathrm{bat}}$'
```

---

## Line Styles and Markers

### Line Specifications
```matlab
% Basic plot with line style
plot(x, y, 'LineSpec', 'LineWidth', 2.5);
```

**Line Styles:**
- `-` : Solid line (default)
- `--` : Dashed line
- `:` : Dotted line
- `-.` : Dash-dot line

**Colors:**
- `'b'` : Blue
- `'r'` : Red
- `'g'` : Green
- `'k'` : Black
- `'m'` : Magenta
- `'c'` : Cyan
- RGB triplet: `[0.2, 0.4, 0.6]`

**Markers:**
- `'o'` : Circle
- `'s'` : Square
- `'^'` : Upward triangle
- `'d'` : Diamond
- `'pentagram'` : Five-pointed star

### Complete Line Specification
```matlab
plot(alpha_deg, charging_time, 'r-', ...
     'LineWidth', 2.5, ...
     'Marker', 's', ...
     'MarkerSize', 8, ...
     'MarkerFaceColor', 'r');
```

### Multiple Line Styles for Comparison
```matlab
hold on;
plot(x, y1, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(x, y2, 'r-s', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(x, y3, 'g-d', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(x, y4, 'k--', 'LineWidth', 3, 'Marker', 'pentagram', 'MarkerSize', 10);
hold off;
```

---

## Axis Configuration

### Font Sizes (Research Paper Standards)
- **Axis labels**: 14-16pt
- **Title**: 18pt
- **Legend**: 12-14pt
- **Tick labels**: 12pt (via `set(gca, 'FontSize', 12)`)

### Setting Axis Properties
```matlab
% Configure current axes
set(gca, 'FontSize', 14, 'LineWidth', 1.2);

% Axis limits
xlim([0, 180]);
ylim([0, max(y)*1.1]);

% Grid
grid on;
```

### Logarithmic Scales
```matlab
% Y-axis logarithmic
semilogy(x, y, 'b-', 'LineWidth', 2);

% X-axis logarithmic
semilogx(x, y, 'b-', 'LineWidth', 2);

% Both axes logarithmic
loglog(x, y, 'b-', 'LineWidth', 2);
```

---

## Legends

### Basic Legend
```matlab
legend('Data 1', 'Data 2', ...
       'Interpreter', 'latex', ...
       'FontSize', 14, ...
       'Location', 'best');
```

### Legend Locations
- `'best'` : Automatic placement
- `'northeast'` : Top right (default)
- `'northwest'` : Top left
- `'southeast'` : Bottom right
- `'southwest'` : Bottom left
- `'north'`, `'south'`, `'east'`, `'west'`

### Legend with LaTeX Math
```matlab
leg_entries = {
    '$P_{\mathrm{batt}} = I_{\mathrm{rms}}^2 R_{\mathrm{bat}}$',
    '$P_{\mathrm{th,cond}} = V_t I_{\mathrm{avg}} + R_{\mathrm{th}} I_{\mathrm{rms}}^2$',
    '$P_{\mathrm{total}}$'
};
legend(leg_entries, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
```

---

## Color Schemes

### Standard Color Palette
For consistent, professional plots:

```matlab
% Primary colors
blue = [0, 0.4470, 0.7410];      % MATLAB default blue
red = [0.8500, 0.3250, 0.0980];  % MATLAB default red
green = [0.4660, 0.6740, 0.1880]; % MATLAB default green
```

### Multi-line Color Scheme
```matlab
colors = lines(7);  % Get 7 distinct colors from 'lines' colormap
plot(x, y1, 'Color', colors(1,:), 'LineWidth', 2);
plot(x, y2, 'Color', colors(2,:), 'LineWidth', 2);
```

---

## Scale Types

### When to Use Different Scales

**Linear Scale** (default):
- Data spans a reasonable range (e.g., 0-100)
- No extreme variations

**Logarithmic Scale** (`semilogy`, `semilogx`, `loglog`):
- Data spans multiple orders of magnitude
- Exponential growth or decay
- Frequency response plots

**Example: Charging Time vs Firing Angle**
```matlab
% Linear scale - good for moderate variations
figure('Position', [100, 100, 800, 600]);
plot(alpha_deg, charging_time_hours, 'r-s', 'LineWidth', 2.5);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Charging Time (hours)', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

% Logarithmic scale - better for large variations
figure('Position', [100, 100, 800, 600]);
semilogy(alpha_deg, charging_time_hours, 'r-s', 'LineWidth', 2.5);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Charging Time (hours, log scale)', 'Interpreter', 'latex', 'FontSize', 16);
grid on;
```

---

## Best Practices

### 1. Line Width
- **Standard plots**: LineWidth = 2-2.5
- **Emphasis lines** (e.g., total): LineWidth = 3
- **Reference lines**: LineWidth = 1.5

### 2. Marker Sizes
- **Standard markers**: MarkerSize = 6-8
- **Important points**: MarkerSize = 10-12

### 3. Font Sizes
```matlab
% Recommended hierarchy
title_size = 18;
label_size = 16;
legend_size = 14;
tick_size = 12;
```

### 4. Grid Lines
```matlab
grid on;
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3); % Subtle grid
```

### 5. Saving Figures
```matlab
% For papers - use vector formats
saveas(gcf, 'figure_name.eps', 'epsc');  % EPS format
saveas(gcf, 'figure_name.pdf', 'pdf');   % PDF format

% For presentations - use high-res raster
print('figure_name.png', '-dpng', '-r300'); % 300 DPI PNG
```

### 6. Aspect Ratio
```matlab
% Maintain aspect ratio
axis equal;  % Equal scaling
axis square; % Square axis box
pbaspect([1.6 1 1]); % Golden ratio
```

---

## Complete Example

```matlab
% Generate data
alpha_deg = 0:5:175;
power_loss = 100 * (1 + cos(deg2rad(alpha_deg)));

% Create figure
figure('Name', 'Power Loss Analysis', 'Position', [100, 100, 900, 600]);

% Plot data
plot(alpha_deg, power_loss, 'b-', ...
     'LineWidth', 2.5, ...
     'Marker', 'o', ...
     'MarkerSize', 8, ...
     'MarkerFaceColor', 'b');

% Configure axes
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 1.2);
xlabel('Firing Angle $\alpha$ (degrees)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Power Loss (W)', 'Interpreter', 'latex', 'FontSize', 16);
title('Power Loss vs Firing Angle', 'Interpreter', 'latex', 'FontSize', 18);

% Add legend
legend('$P_{\mathrm{loss}}$', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% Set axis limits
xlim([0, 180]);
ylim([0, max(power_loss)*1.1]);

% Save figure
saveas(gcf, 'power_loss.pdf', 'pdf');
```

---

## References

- [MATLAB Graphics Documentation](https://www.mathworks.com/help/matlab/graphics.html)
- [LaTeX Mathematical Symbols](https://www.overleaf.com/learn/latex/List_of_Greek_letters_and_math_symbols)
- [MATLAB Plot Customization](https://www.mathworks.com/help/matlab/creating_plots/customize-graph-appearance.html)

---

**Version**: 1.0  
**Last Updated**: November 2025  
**Project**: ThyristoLab - Battery Charger Analysis
