% Create right curve

[probability_estimate,xi] = ksdensity((signal)); % get the densities

figure(2)
plot(xi,probability_estimate,'-k','LineWidth',2);
ylabel('Density'); xlabel('log2(FPKM)')
hold on


peak_max_idx = (find(probability_estimate==max(probability_estimate))); % find the maximum of the density curve
max_r_side = probability_estimate(peak_max_idx:end); % copy right-most side of density curve
max_l_side = flip(max_r_side); % mirror the right side
hybrid = [max_l_side(1:end-1), max_r_side];

hybrid_curve = zeros(numel(probability_estimate),1); % create hybrid curve
if numel(hybrid)> numel(probability_estimate)
    hybrid_curve(end+1-numel(hybrid(end-100+1:end)):end)=hybrid(end-100+1:end);
    x = zeros(numel(probability_estimate),1); % create new x-axis
else
    hybrid_curve(end-numel(hybrid)+1:end) = hybrid;
    x = xi;
end


plot(x,hybrid_curve,'b--','LineWidth',2);
legend({'Signal','hybrid curve right'});
title({'Density plot with hybrid curves',['Sample: ', colnames{1,1}]},'interpreter','none');

clear peak_max_idx max_l_side max_r_side hybrid
%%
% Create left curve (main curve - rightmost curve)

rest_curve = probability_estimate - hybrid_curve';
rest_curve(find(rest_curve<0.0001))=0;


plot(xi,rest_curve, 'r--','LineWidth',2);
legend({'Signal','hybrid curve right','hybrid curve left'});









%% fit the curves

