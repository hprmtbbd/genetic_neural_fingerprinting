function [cmat, tpr, tnr, ppv, AUC] = calc_ML_stats(model,true_labels,pred_labels,pred_scores)
cmat = confusionmat(true_labels,pred_labels,'Order',model.ClassNames);

n = []; tpr = []; tnr = []; ppv = []; AUC = [];
for j = 1:size(pred_scores,2)
    n(j) = sum(cmat(j,:));
    tpr(j) = cmat(j,j)/sum(cmat(j,:));
    
    cmat_temp = cmat; cmat_temp(j,:) = [];
    tnr(j) = 1-(sum(cmat_temp(:,j))/sum(sum(cmat_temp)));
    
    ppv(j) = cmat(j,j)/sum(cmat(:,j));
    
    [~,~,~,AUC(j)] = perfcurve(true_labels,pred_scores(:,j),model.ClassNames{j});
end

wmean = @(x,n) (x*n')/sum(n);
tpr = [tpr, wmean(tpr,n)]; tnr = [tnr, wmean(tnr,n)]; ppv = [ppv, wmean(ppv,n)]; AUC = [AUC, wmean(AUC,n)];

end