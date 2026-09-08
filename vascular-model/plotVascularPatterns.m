function plotVascularPatterns(t, y, g, params, tIdx)
% PLOTVASCULARPATTERNS  Visualize the Carteni et al. vascular pattern simulation.
%
%   plotVascularPatterns(t, y, g, params)
%       Plots final state of all species.
%
%   plotVascularPatterns(t, y, g, params, tIdx)
%       Plots state at specified time indices.

    if nargin < 5
        tIdx = length(t);  % default: final state
    end

    N = params.gridSize;
    nTot = N * N;
    mask = g.mask;

    % Offsets for each species in the state vector
    off_S0 = 0;  off_S1 = 1;  off_S2 = 2;
    off_AP = 3;  off_AF = 4;  off_AX = 5;  off_R  = 6;

    % Helper to extract and reshape a species field
    function F = getField(yRow, offset)
        F = reshape(yRow(offset*nTot + (1:nTot)), [N, N]);
    end

    % Build figure
    nTimes = numel(tIdx);
    nRows = nTimes;
    nCols = 7;

    figure('Name', 'Vascular Pattern Simulation', ...
           'Position', [50, 50, 1800, 250 * nRows], ...
           'Color', 'w');

    speciesNames = {'S_0 (morphogen)', 'S_1 (central subst.)', 'S_2 (periph. subst.)', ...
                    'A_P (procambium)', 'A_F (phloem)', 'A_X (xylem)', 'R (repressor)'};
    offsets      = [off_S0, off_S1, off_S2, off_AP, off_AF, off_AX, off_R];
    cmaps        = {@parula, @hot, @cool, @summer, @autumn, @winter, @gray};

    for row = 1:nTimes
        ti = tIdx(row);
        yRow = y(ti, :);

        for col = 1:7
            ax = subplot(nRows, nCols, (row-1)*nCols + col);
            F = getField(yRow, offsets(col));
            F(~mask) = NaN;

            imagesc(g.X(1,:), g.Y(:,1), F);
            axis equal tight;
            colormap(ax, cmaps{col}(256));
            colorbar;
            set(gca, 'YDir', 'normal');

            if row == 1
                title(speciesNames{col}, 'FontSize', 12, 'FontWeight', 'bold');
            end

            if col == 1
                ylabel(sprintf('t = %.1f', t(ti)), 'FontSize', 11);
            end

            if row == nRows
                xlabel('x');
                if col == 1, ylabel('y'); end
            end
        end
    end

    sgtitle('Carteni et al. (2014) - Primary Vascular Structure Model', ...
            'FontSize', 14, 'FontWeight', 'bold');

    % ---- Summary figure ----
    figure('Name', 'Vascular Differentiation Summary', ...
           'Position', [100, 100, 1200, 500], 'Color', 'w');

    ti = tIdx(end);
    yRow = y(ti, :);
    AP = getField(yRow, off_AP);
    AF = getField(yRow, off_AF);
    AX = getField(yRow, off_AX);

    % Procambium
    subplot(1, 4, 1);
    F = AP; F(~mask) = NaN;
    imagesc(g.X(1,:), g.Y(:,1), F);
    axis equal tight; colorbar; set(gca, 'YDir', 'normal');
    title('Procambium (A_P)', 'FontSize', 13);
    xlabel('x'); ylabel('y');

    % Phloem
    subplot(1, 4, 2);
    F = AF; F(~mask) = NaN;
    imagesc(g.X(1,:), g.Y(:,1), F);
    axis equal tight; colorbar; set(gca, 'YDir', 'normal');
    title('Phloem (A_F)', 'FontSize', 13);
    xlabel('x'); ylabel('y');

    % Xylem
    subplot(1, 4, 3);
    F = AX; F(~mask) = NaN;
    imagesc(g.X(1,:), g.Y(:,1), F);
    axis equal tight; colorbar; set(gca, 'YDir', 'normal');
    title('Xylem (A_X)', 'FontSize', 13);
    xlabel('x'); ylabel('y');

    % Overlay
    subplot(1, 4, 4);
    overlay = zeros(N, N, 3);
    AF_norm = AF ./ max(AF(mask), [], 'all');
    AX_norm = AX ./ max(AX(mask), [], 'all');
    overlay(:,:,1) = AX_norm;
    overlay(:,:,2) = AF_norm;
    overlay(:,:,3) = 0;
    overlay(repmat(~mask, [1,1,3])) = 1;
    image(g.X(1,:), g.Y(:,1), overlay);
    axis equal tight; set(gca, 'YDir', 'normal');
    title('Overlay: Xylem (red) + Phloem (green)', 'FontSize', 13);
    xlabel('x'); ylabel('y');

    sgtitle(sprintf('Final State (t = %.1f)', t(ti)), ...
            'FontSize', 14, 'FontWeight', 'bold');

    % ---- Radial profiles ----
    figure('Name', 'Radial Profiles', ...
           'Position', [150, 150, 1000, 400], 'Color', 'w');

    centreRow = round(N/2);
    r_vals = g.X(centreRow, :);

    subplot(1, 2, 1);
    plot(r_vals, AP(centreRow, :), 'b-', 'LineWidth', 2); hold on;
    plot(r_vals, AF(centreRow, :), 'g-', 'LineWidth', 2);
    plot(r_vals, AX(centreRow, :), 'r-', 'LineWidth', 2);
    S0 = getField(yRow, off_S0);
    plot(r_vals, S0(centreRow, :), 'k--', 'LineWidth', 1.5);
    xline(0, 'k:'); yline(params.AP_star, 'b:');
    xlabel('x position'); ylabel('Concentration');
    legend('A_P (procambium)', 'A_F (phloem)', 'A_X (xylem)', 'S_0 (morphogen)', ...
           'Location', 'best');
    title('Horizontal Cross-Section Profiles', 'FontSize', 12);
    grid('on');

    % Radial scatter
    subplot(1, 2, 2);
    [sortedD, sortIdx] = sort(g.D(mask));
    AP_vec = AP(mask); AF_vec = AF(mask); AX_vec = AX(mask); S0_vec = S0(mask);
    plot(sortedD, AP_vec(sortIdx), 'b.', 'MarkerSize', 4); hold on;
    plot(sortedD, AF_vec(sortIdx), 'g.', 'MarkerSize', 4);
    plot(sortedD, AX_vec(sortIdx), 'r.', 'MarkerSize', 4);
    plot(sortedD, S0_vec(sortIdx), 'k.', 'MarkerSize', 2);
    xlabel('Distance from centre'); ylabel('Concentration');
    legend('A_P', 'A_F', 'A_X', 'S_0', 'Location', 'best');
    title('Radial Distribution', 'FontSize', 12);
    grid('on');

    drawnow;
end
