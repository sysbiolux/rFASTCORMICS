function [discretized] = discretize_rpkm(rpkm,colnames,figflag)
%% Discretize_Data
% script adapted from (c) Dr. Maria Pires Pacheco 2016

%% log2-transform the data
signal = log2(rpkm);
signal(isinf(signal)) = -10000;

% figflag = 0; %show figures
%% Discretize the data by creating half-gaussians
for j=1:size(rpkm,2) %for each sample
    signal_sample = signal(:,j);
    data_keep = rpkm(:,j);
        signal_sample = signal_sample(exp(signal_sample) > 0);

%     signal_sample = signal_sample(signal_sample > -max(signal_sample));
%     signal_sample = signal_sample(signal_sample > -100);
    
    % Density plot
    [probability_estimate,xi] = ksdensity((signal_sample));
    if figflag;figure; hold on
        plot(xi,probability_estimate,'-k','LineWidth',2);
        title({'Density curve for sample: ', colnames{j}},'Interpreter','none');
        ylabel('Density'); xlabel('log2(rpkm)');
        title({'Density plot',['Sample: ', colnames{j}]},'Interpreter','none');
        legend({'Signal'},'Location','best')
%         saveas(gcf,['Figures\density_', colnames{j},'.png']);
    end
    %% create right curve
    
    peak_max_idx = (find(probability_estimate==max(probability_estimate))); % find the maximum of the density curve
    max_r_side = probability_estimate(peak_max_idx:end); % copy right-most side of density curve
    max_l_side = flip(max_r_side); % mirror the right side
    hybrid = [max_l_side(1:end-1), max_r_side];
    
    
    hybrid_curve = zeros(numel(probability_estimate),1);
    if numel(hybrid)> numel(probability_estimate)
        hybrid_curve(end+1-numel(hybrid(end-100+1:end)):end)=hybrid(end-100+1:end);
%         x = zeros(numel(probability_estimate),1); % create mew x-axis
        x = xi;

    else
        hybrid_curve(end-numel(hybrid)+1:end) = hybrid;
        x = xi;
    end
    
    if figflag;
        plot(x,hybrid_curve,'b--','LineWidth',2);
        legend({'Signal','hybrid curve right'},'Location','best')
        title({'Density plot with hybrid curves',['Sample: ', colnames{j}]},'Interpreter','none');end;
    
    clear peak_max_idx max_l_side max_r_side hybrid
    %% create left curve (main curve - rightmost curve)
    
    rest_curve = probability_estimate - hybrid_curve';
    rest_curve(find(rest_curve<0.0001))=0;
    
    if figflag; plot(xi,rest_curve, 'r--','LineWidth',2);
        legend({'Signal','hybrid curve right','hybrid curve left'},'Location','best')
        saveas(gcf,['Figures\density_maxfit_', colnames{j},'.png']);
    end
    %% fit the curves
    
    [fitresult_r, ~] = createFit_r(xi, hybrid_curve);
    [fitresult_l, a] = createFit_l(xi, rest_curve);
    
    if figflag;
        plot(fitresult_r,'b');
        plot(fitresult_l),'g';
        title({'Density plot with fitted curves',['Sample: ', colnames{j}]},'Interpreter','none');
        ylabel('Density'); xlabel('log2(rpkm)');
        legend({'Signal', 'hybrid curve right', 'hybrid curve left','fitted curve right','fitted curve left'},'Location','best')
%         saveas(gcf,['Figures\density_fitted_', colnames{j},'.png']);
    end
    %% zFRMA transform the data and plot the density
    
    sigma1 = fitresult_r.c1/sqrt(2);
    mu1 = fitresult_r.b1; % x-value of right max
    
    zFRMA = (signal_sample-mu1)/sigma1;
    [yFRMA,xFRMA] = ksdensity(zFRMA);
    
    
    clear hybrid_curve rest_curve yFRMA xFRMA x xi ans
    %% do the analysis
    
    zFRMA= (signal(:,j)-mu1)/sigma1;
    
    sigma2 = fitresult_l.c1/sqrt(2);
    mu2 = fitresult_l.b1; % x-value of left max
    mu2_t=(mu2-mu1)/sigma1;
    
    if figflag
        plot([mu1, mu1],[0, max(probability_estimate)])
        plot([mu2, mu2],[0, max(probability_estimate)])
        legend({'Signal', 'hybrid curve right', 'hybrid curve left',...
            'fitted curve right','fitted curve left',...
            ['expression threshold = ',num2str(mu1)],['inexpression threshold = ',num2str(mu2)]},'Location','best')
        saveas(gcf,['Figures\density_fitted_mu_', colnames{j},'.png']);close all;

    end
    
    discretized=zFRMA;
    
    e=0;
    ue=max (-3, mu2_t);
    zFRMA=reshape(zFRMA,size(data_keep,1),size(data_keep,2));
    
    exp_threshold=e;
    unexp_threshold=ue;
    exp=find(discretized>=exp_threshold);
    not_exp=find(discretized<=unexp_threshold);
    
    unknown=find(discretized<exp_threshold & discretized>unexp_threshold);
    discretized(unknown)=0;
    discretized(not_exp)=-1;
    discretized(exp)=1;
    
    discretized=(reshape(discretized,size(data_keep,1),size(data_keep,2)));
    clear fitresult_l fitresult_r data_keep e exp exp_threshold mu1 mu2 mu2_t not_exp sigma1 sigma2 signal2 signal_sample ue unexp_threshold unknown zFRMA
    discretized_keep(j,:)=discretized';
end
discretized = discretized_keep';
clear discretized_keep j
clear j probability_estimate signal_sample xi