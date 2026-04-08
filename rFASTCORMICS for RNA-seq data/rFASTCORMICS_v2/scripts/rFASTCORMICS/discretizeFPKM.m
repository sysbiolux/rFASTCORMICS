function [discretized] = discretizeFPKM(fpkm, colnames, figflag, pathFigures)
%   The function discretizes gene expression data (FPKM values) into three categories:
%   expressed, not expressed, and unknown, based on a zFPKM transformation and
%   half-Gaussian density fitting.
%
% USAGE:
%
%   [discretized] = discretizeFPKM(fpkm, colnames, figflag, pathFigures)
%
% INPUTS:
%   fpkm:              m x n matrix or table of FPKM expression values (genes x samples)
%   colnames:          1 x n cell array of sample names
%   figflag:           optional, flag to plot and save figures (default: 0)
%   pathFigures:       optional, path to save figures (default: current folder)
%
% OUTPUT:
%   discretized:       m x n matrix of discretized values
%                      *  1 : expressed
%                      *  0 : unknown
%                      * -1 : not expressed
%
% .. Authors:
%       - Maria Pires Pacheco, 2016, University of Luxembourg

if nargin < 3
    figflag = 0;
end
if nargin < 4
    pathFigures = '';
end

%% Convert table input to array if needed
if istable(fpkm)
    fpkm = table2array(fpkm);
end

%% log2-transform the data
signal = log2(fpkm);
signal(isinf(signal)) = -1e6;

%% Discretize the data by creating half-gaussians
for j = 1:size(fpkm, 2) % for each sample
    signal_sample = signal(:, j);
    data_keep = fpkm(:, j);
    signal_sample = signal_sample(signal_sample > -1e6);
    
    % Histogram
    if figflag
        figure
        % create folder to save Figures if doesnt exist
        if ~strcmp(pathFigures, '')
            if ~exist([pathFigures, 'Figures/Discretization/'], 'dir')
                mkdir([pathFigures, 'Figures/Discretization/'])
            end
        else
            if ~exist('Figures/Discretization/', 'dir')
                mkdir('Figures/Discretization/')
            end
        end
        
        histogram((signal_sample), 100);
        ylabel('abundance');
        xlabel('log2(FPKM)');
        title({'Histogram of abundance', ['Sample: ', colnames{j}]}, 'Interpreter', 'none');
        if strcmp(pathFigures, '')
            saveas(gcf, ['Figures\Discretization\histogram_', colnames{j}, '.png'])
        else
            saveas(gcf,[pathFigures,'Figures\Discretization\histogram_', colnames{j},'.png'])
        end
    end
    
    % Density plot
    [probability_estimate, xi] = ksdensity((signal_sample));
    if figflag
        figure; 
        hold on;
        plot(xi, probability_estimate, '-k', 'LineWidth', 2);
        ylabel('Density'); xlabel('log2(FPKM)');
        title({'Density plot', ['Sample: ', colnames{j}]}, 'Interpreter', 'none');
        legend({'Signal'}, 'Location', 'best')
        saveas(gcf, [pathFigures, 'Figures\Discretization\density_', colnames{j}, '.png']);
    end

    %% create right curve
    peak_max_idx = (find(probability_estimate == max(probability_estimate))); % find the maximum of the density curve
    max_r_side = probability_estimate(peak_max_idx:end); % copy right-most side of density curve
    max_l_side = flip(max_r_side); % mirror the right side
    hybrid = [max_l_side(1:end-1), max_r_side];
    
    % Adjust hybrid curve to match length of density
    hybrid_curve = zeros(numel(probability_estimate), 1);
    if numel(hybrid) > numel(probability_estimate)
        hybrid_curve(end+1-numel(hybrid(end-100+1:end)):end) = hybrid(end-100+1:end);
        % x = zeros(numel(probability_estimate),1); % create new x-axis
        x = xi;
    else
        hybrid_curve(end-numel(hybrid)+1:end) = hybrid;
        x = xi;
    end
    
    if figflag
        plot(x, hybrid_curve, 'b--', 'LineWidth', 2);
        legend({'Signal','hybrid curve right'}, 'Location', 'best')
        title({'Density plot with hybrid curves', ['Sample: ', colnames{j}]}, 'Interpreter', 'none');
    end
    
    %% create left curve (main curve - rightmost curve)
    rest_curve = probability_estimate - hybrid_curve';
    rest_curve(rest_curve < 0.0001) = 0;
    
    if figflag
        plot(xi, rest_curve, 'r--', 'LineWidth', 2);
        legend({'Signal', 'hybrid curve right', 'hybrid curve left'}, 'Location', 'best')
        saveas(gcf, [pathFigures,'Figures\Discretization\density_maxfit_', colnames{j}, '.png']);
    end

    %% fit half-Gaussian curves to right and left curves
    [fitresult_r, ~] = createFitR(xi, hybrid_curve);
    [fitresult_l, ~] = createFitL(xi, rest_curve);
    
    if figflag
        plot(fitresult_r, 'b');
        plot(fitresult_l, 'g');
        title({'Density plot with fitted curves', ['Sample: ', colnames{j}]}, 'Interpreter','none');
        ylabel('Density');
        xlabel('log2(FPKM)');
        legend({'Signal', 'hybrid curve right', 'hybrid curve left', 'fitted curve right', 'fitted curve left'}, 'Location','best')
        saveas(gcf, [pathFigures, 'Figures\Discretization\density_fitted_', colnames{j}, '.png']);
        close all;
    end
    
    %% zFPKM transform the data and plot the density
    
    sigma1 = fitresult_r.c1/sqrt(2);
    mu1 = fitresult_r.b1; % x-value of right max
    
    zFPKM = (signal_sample - mu1)/sigma1;
    [yFPKM, xFPKM] = ksdensity(zFPKM);
    
    %% do the analysis
    
    zFPKM = (signal(:, j) - mu1)/sigma1; %transform sample into zFPKM
    
    mu2 = fitresult_l.b1;
    mu2_t = (mu2-mu1)/sigma1; % x-value of left max
    
    if figflag
        plot([mu1, mu1], [0, max(probability_estimate)])
        plot([mu2, mu2], [0, max(probability_estimate)])
        legend({'Signal', 'hybrid curve right', 'hybrid curve left',...
            'fitted curve right', 'fitted curve left',...
            ['expression threshold = ', num2str(mu1)], ['inexpression threshold = ', num2str(mu2)]}, 'Location', 'best')
    end
    
    discretized = zFPKM;
    
    e = 0;
    ue = max(-3, mu2_t);
    
    if figflag
        figure;
        hold on;
        plot(xFPKM, yFPKM, '-k', 'LineWidth', 2);
        xlabel('zFPKM');
        ylabel('density');
        line([0 0], [0.0 0.3],'color', [0, 1, 0])
        line([ue ue], [0.0 0.3],'color', [0, 1, 1]);
        legend({'zFPKM data', 'expression threshold', 'inexpression threshold'}, 'Location', 'best')
        title({'Expression thresholds', ['Sample: ', colnames{j}], ['inexpression threshold = ', num2str(ue)]}, 'Interpreter', 'none');
        
        saveas(gcf, [pathFigures, 'Figures\Discretization\zFPKM_', colnames{j}, '.png']);
        close all;
    end
    
    exp_threshold = e;
    unexp_threshold = ue;
    
    expressed = find(discretized >= exp_threshold);
    not_exp = find(discretized <= unexp_threshold);
    unknown = find(discretized < exp_threshold & discretized > unexp_threshold);
    
    discretized(unknown) = 0;
    discretized(not_exp) = -1;
    discretized(expressed) = 1;
    
    discretized = (reshape(discretized, size(data_keep, 1), size(data_keep, 2)));
    discretized_keep(j, :) = discretized';
end
% Combine all samples
discretized = discretized_keep';

%% Density plot for all samples
if figflag
    figure; 
    hold on;
    for j = 1:length(colnames) %for each sample
        signal_sample = signal(:, j);
        signal_sample = signal_sample(exp(signal_sample) > 1e-6);
        
        % Density plot
        [probability_estimate, xi] = ksdensity((signal_sample));
        % plot(xi, probability_estimate, line_col{j}, 'LineWidth', 1); end
        plot(xi, probability_estimate, 'LineWidth', 1);
    end
    
    title({'Density curves for all samples'});
    ylabel('Density'); 
    xlabel('log2(FPKM)');
    % legend(colnames,'interpreter','none','Location','best')
    if strcmp(pathFigures, '')
        saveas(gcf, 'Figures\Discretization\density_all.png'); 
        close all
    else
        saveas(gcf,[pathFigures,'Figures\Discretization\density_all.png']); 
        close all
    end
end

%% Subplots for all samples
if figflag
    figure; 
    hold on;
    for j = 1:length(colnames) %for each sample
        signal_sample = signal(:, j);
        signal_sample = signal_sample(exp(signal_sample) > 1e-6);
        
        % Density plot
        [probability_estimate, xi] = ksdensity((signal_sample));
        subplot(5, 6, j)
        % plot(xi, probability_estimate, line_col{j}, 'LineWidth', 1);
        plot(xi, probability_estimate, 'LineWidth', 1);
        title({colnames{j}}, 'interpreter', 'none');
        % ylabel('Density'); 
        % xlabel('log2(FPKM)');
    end
    if strcmp(pathFigures, '')
        saveas(gcf, 'Figures\Discretization\density_plots.png'); 
        close all
    else
        saveas(gcf, [pathFigures, 'Figures\Discretization\density_plots.png']); 
        close all
    end
end
