function [cvm_crit,acc,cmat,tpr,tnr,ppv,auc,sel_featnames] = run_fingerprint_train_test(dat,covar_dat,iall,jall,train_ind,test_ind,all_label,feat_names,analysis_flag,measure_flag,covar_flag,fs_flag,diff_func0,dist_func0,fitmodel,input_add)

%%%%
jtrain = jall(train_ind);
itrain = iall(train_ind);
train_label = all_label(train_ind);

jtrain_ur = jtrain(strcmp(train_label,'UR'));
itrain_ur = itrain(strcmp(train_label,'UR'));

jtrain_mz = jtrain(strcmp(train_label,'MZ'));
itrain_mz = itrain(strcmp(train_label,'MZ'));

jtrain_sib = jtrain(strcmp(train_label,'SIB'));
itrain_sib = itrain(strcmp(train_label,'SIB'));

if analysis_flag == 1
    dat_train_ind = [jtrain_mz; itrain_mz; jtrain_sib; itrain_sib; jtrain_ur; itrain_ur];
elseif analysis_flag == 2
    dat_train_ind = [jtrain_mz; itrain_mz; jtrain_ur; itrain_ur];
elseif analysis_flag == 3
    dat_train_ind = [jtrain_mz; itrain_mz; jtrain_sib; itrain_sib];
elseif analysis_flag == 4
    dat_train_ind = [jtrain_sib; itrain_sib; jtrain_ur; itrain_ur];
end
u_dat_train_ind = unique(dat_train_ind,'stable');

if measure_flag == 1
    dat_train = dat(u_dat_train_ind,:);
    
    if covar_flag
        covar_train = covar_dat(u_dat_train_ind,:);
        Dsn_model = [0 0; 1 0; 0 1; 0 2; 1 1; 1 2];
        design_mat = x2fx(covar_train,Dsn_model,1);
        beta = design_mat\dat_train;
        dat_train_out = dat_train - design_mat*beta;
    else
        dat_train_out = dat_train;
    end
    
    mdat = mean(dat_train_out);
    sigma = std(dat_train_out);
    dat_train0 = (dat_train_out - mdat)./sigma;
    
    pdat_train = dat_train0;
    
elseif measure_flag == 2
    dat_train = dat(u_dat_train_ind,:,:);
    dat_train_out = dat_train;
    
    mdat = mean(dat_train_out);
    sigma = std(dat_train_out);
    dat_train0 = (dat_train_out - mdat)./sigma;
    
    pdat_train = dat_train0;
end

if ismember(analysis_flag,[1,2,3])
    [~,jmz0] = ismember(jtrain_mz,u_dat_train_ind);
    [~,imz0] = ismember(itrain_mz,u_dat_train_ind);
    if measure_flag == 1
        dtrain_mz = diff_func0(pdat_train(jmz0,:),pdat_train(imz0,:));
    elseif measure_flag == 2
        dtrain_mz = dist_func0(pdat_train(jmz0,:,:),pdat_train(imz0,:,:));
    end
end

if ismember(analysis_flag,[1,3,4])
    [~,jsib0] = ismember(jtrain_sib,u_dat_train_ind);
    [~,isib0] = ismember(itrain_sib,u_dat_train_ind);
    if measure_flag == 1
        dtrain_sib = diff_func0(pdat_train(jsib0,:),pdat_train(isib0,:));
    elseif measure_flag == 2
        dtrain_sib = dist_func0(pdat_train(jsib0,:,:),pdat_train(isib0,:,:));
    end
end

if ismember(analysis_flag,[1,2,4])
    [~,jur0] = ismember(jtrain_ur,u_dat_train_ind);
    [~,iur0] = ismember(itrain_ur,u_dat_train_ind);
    if measure_flag == 1
        dtrain_ur = diff_func0(pdat_train(jur0,:),pdat_train(iur0,:));
    elseif measure_flag == 2
        dtrain_ur = dist_func0(pdat_train(jur0,:,:),pdat_train(iur0,:,:));
    end
end

if analysis_flag == 1
    dtrain = [dtrain_mz; dtrain_sib; dtrain_ur];
elseif analysis_flag == 2
    dtrain = [dtrain_mz; dtrain_ur];
elseif analysis_flag == 3
    dtrain = [dtrain_mz; dtrain_sib];
elseif analysis_flag == 4
    dtrain = [dtrain_sib; dtrain_ur];
end

if analysis_flag == 1
    train_label_t = train_label;
elseif analysis_flag == 2
    train_label_t = train_label(~strcmp(train_label,'SIB'));
elseif analysis_flag == 3
    train_label_t = train_label(~strcmp(train_label,'UR'));
elseif analysis_flag == 4
    train_label_t = train_label(~strcmp(train_label,'MZ'));
end

%% featsel
if fs_flag == 0
    cvm_crit = NaN;
    sel_feat = 1:size(dtrain,2);
end
sel_featnames = feat_names(sel_feat);

dtrain_in = dtrain(:,sel_feat);

input = [{dtrain_in,train_label_t},input_add];
model = fitmodel(input{:});

%%
jtest = jall(test_ind);
itest = iall(test_ind);
test_label = all_label(test_ind);

jtest_ur = jtest(strcmp(test_label,'UR'));
itest_ur = itest(strcmp(test_label,'UR'));

jtest_mz = jtest(strcmp(test_label,'MZ'));
itest_mz = itest(strcmp(test_label,'MZ'));

jtest_sib = jtest(strcmp(test_label,'SIB'));
itest_sib = itest(strcmp(test_label,'SIB'));

if analysis_flag == 1
    dat_test_ind = [jtest_mz; itest_mz; jtest_sib; itest_sib; jtest_ur; itest_ur];
elseif analysis_flag == 2
    dat_test_ind = [jtest_mz; itest_mz; jtest_ur; itest_ur];
elseif analysis_flag == 3
    dat_test_ind = [jtest_mz; itest_mz; jtest_sib; itest_sib];
elseif analysis_flag == 4
    dat_test_ind = [jtest_sib; itest_sib; jtest_ur; itest_ur];
end
u_dat_test_ind = unique(dat_test_ind,'stable');

if measure_flag == 1
    dat_test = dat(u_dat_test_ind,:);
    
    if covar_flag
        covar_test = covar_dat(u_dat_test_ind,:);
        Dsn_model = [0 0; 1 0; 0 1; 0 2; 1 1; 1 2];
        design_mat = x2fx(covar_test,Dsn_model,1);
        dat_test_out = dat_test - design_mat*beta;
    else
        dat_test_out = dat_test;
    end
    
    dat_test0 = (dat_test_out - mdat)./sigma;
    pdat_test = dat_test0;
elseif measure_flag == 2
    dat_test = dat(u_dat_test_ind,:,:);
    dat_test_out = dat_test;
    
    dat_test0 = (dat_test_out - mdat)./sigma;
    pdat_test = dat_test0;
end

if ismember(analysis_flag,[1,2,3])
    [~,jmz0] = ismember(jtest_mz,u_dat_test_ind);
    [~,imz0] = ismember(itest_mz,u_dat_test_ind);
    if measure_flag == 1
        dtest_mz = diff_func0(pdat_test(jmz0,:),pdat_test(imz0,:));
    elseif measure_flag == 2
        dtest_mz = dist_func0(pdat_test(jmz0,:,:),pdat_test(imz0,:,:));
    end
end

if ismember(analysis_flag,[1,3,4])
    [~,jsib0] = ismember(jtest_sib,u_dat_test_ind);
    [~,isib0] = ismember(itest_sib,u_dat_test_ind);
    if measure_flag == 1
        dtest_sib = diff_func0(pdat_test(jsib0,:),pdat_test(isib0,:));
    elseif measure_flag == 2
        dtest_sib = dist_func0(pdat_test(jsib0,:,:),pdat_test(isib0,:,:));
    end
end

if ismember(analysis_flag,[1,2,4])
    [~,jur0] = ismember(jtest_ur,u_dat_test_ind);
    [~,iur0] = ismember(itest_ur,u_dat_test_ind);
    if measure_flag == 1
        dtest_ur = diff_func0(pdat_test(jur0,:),pdat_test(iur0,:));
    elseif measure_flag == 2
        dtest_ur = dist_func0(pdat_test(jur0,:,:),pdat_test(iur0,:,:));
    end
end

if analysis_flag == 1
    dtest = [dtest_mz; dtest_sib; dtest_ur];
elseif analysis_flag == 2
    dtest = [dtest_mz; dtest_ur];
elseif analysis_flag == 3
    dtest = [dtest_mz; dtest_sib];
elseif analysis_flag == 4
    dtest = [dtest_sib; dtest_ur];
end

if analysis_flag == 1
    test_label_t = test_label;
elseif analysis_flag == 2
    test_label_t = test_label(~strcmp(test_label,'SIB'));
elseif analysis_flag == 3
    test_label_t = test_label(~strcmp(test_label,'UR'));
elseif analysis_flag == 4
    test_label_t = test_label(~strcmp(test_label,'MZ'));
end

dtest_in = dtest(:,sel_feat);

%% testing
acc = (1-loss(model,dtest_in,test_label_t));
[pred_labels,pred_score] = predict(model,dtest_in);
[cmat, tpr, tnr, ppv, auc] = calc_ML_stats(model,test_label_t,pred_labels,pred_score);

end

function model1 = fitmodel1(Xtrain,Ytrain,fitmodel,input_add)

input = [{Xtrain,Ytrain},input_add];
model1 = fitmodel(input{:});
    
end
