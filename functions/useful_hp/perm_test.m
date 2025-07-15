function [sig,px_a,pj,xj] = perm_test(metric_true,metric_perm,alpha,side_flag)

niter = length(metric_perm);
[p_dist,x_dist] = ecdf(metric_perm);
% plot(x_dist,p_dist)

% https://www.mathworks.com/help/stats/examples/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html
ndis = length(x_dist)-1;
xj = x_dist(2:end);
pj = (p_dist(1:end-1)+p_dist(2:end))/2;
xj = [xj(1)-pj(1)*(xj(2)-xj(1))/((pj(2)-pj(1))); xj; xj(ndis)+(1-pj(ndis))*((xj(ndis)-xj(ndis-1))/(pj(ndis)-pj(ndis-1)))];
pj = [0; pj; 1];

try
    px = interp1(xj,pj,metric_true,'linear','extrap'); % left-sided p-value
catch e
    % F = ksdensity(sort(metric_perm),sort(metric_perm),'function','cdf');
    % plot(sort(metric_perm),F),hold on,plot(x_dist,p_dist)
    px = ksdensity(sort(metric_perm),metric_true,'function','cdf');
end

px_a = px;
if px_a < 1/niter
    px_a = 1/niter;
elseif px_a > 1-(1/niter)
    px_a = 1-(1/niter);
end

if side_flag == 1
    % right-sided test
    px_a = 1-px_a;
    
    hloc = find(p_dist >= (1-alpha) - eps);
    hval = x_dist(hloc(1));
    if metric_true > hval, sig = 1; else, sig = 0; end
elseif side_flag == -1
    % left-sided test
    
    lloc = find(p_dist <= alpha + eps);
    lval = x_dist(lloc(end));
    if metric_true < lval, sig = 1; else, sig = 0; end
elseif side_flag == 0
    % two-sided test
    if px_a > 0.5, px_a = 1-px_a; end
    px_a = 2*px_a;
    
    balpha = alpha/2;
    hloc = find(p_dist >= (1-balpha) - eps);
    hval = x_dist(hloc(1));
    lloc = find(p_dist <= balpha + eps);
    lval = x_dist(lloc(end));
    if metric_true > hval || metric_true < lval, sig = 1; else, sig = 0; end
end

end
