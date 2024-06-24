function ML_STATS = run_fingerprint(dat,covar_dat,feat_names,imz,jmz,isib,jsib,phi,fam_locs,fam_locs00,nsubj,nmz,nsib,analysis_flag,data_flag,measure_flag,covar_flag,nfold,niter,fs_flag,diff_func0,dist_func0,fitmodel,input_add)

cvm_crit_all = [];
acc_all = [];
cmat_all = []; tpr_all = []; tnr_all = []; ppv_all = []; auc_all = [];
sel_feats_all = {};

parfor itr = 1:niter % parfor
    disp([num2str(itr),'/',num2str(niter)])
    
    if data_flag == 1
        % randomly pair unrelated subjects
        
        nur = round((nmz+nsib)/2);
        
        jur = randsample((1:nsubj)',nur);
        iur = [];
        for i_jur = 1:length(jur)
            jur_i = jur(i_jur);
            phi0_b_vec = phi(:,jur_i);
            iur_ii = find(phi0_b_vec==0);
            iur_i = randsample(iur_ii,1);
            iur(i_jur,1) = iur_i;
        end
        
        npair = nmz+nsib+nur;
        
        jall = [jmz; jsib; jur];
        iall = [imz; isib; iur];
        
        all_label = cell(npair,1);
        all_label(1:nmz) = {'MZ'};
        all_label(nmz+1:nmz+nsib) = {'SIB'};
        all_label(nmz+nsib+1:npair) = {'UR'};
        
    elseif data_flag == 2
        % 1) randomly pair unrelated subjects
        % 2) randomly select sibling pairs from non-MZ families
        % 3) randomly sample from majority class (sibling pairs)
        
        % 1)
        nur = nmz;
        
        jur = randsample((1:nsubj)',nur);
        iur = [];
        for i_jur = 1:length(jur)
            jur_i = jur(i_jur);
            fam_lia = cell2mat(cellfun(@(x) ismember(jur_i,x),fam_locs,'UniformOutput',0));
            not_fam_vec = fam_locs(~fam_lia)';
            not_fam_inds = [not_fam_vec{:}]';
            iur_i = randsample(not_fam_inds,1);
            iur(i_jur,1) = iur_i;
        end
        
        % 2)
        isib_all = []; jsib_all = [];
        for ifam = 1:length(fam_locs00)
            fam_loc = fam_locs00{ifam};
            perm_ind0 = randsample(length(fam_loc),2);
            fam_loc_perm = fam_loc(perm_ind0);
            isib_all(ifam,1) = fam_loc_perm(1);
            jsib_all(ifam,1) = fam_loc_perm(2);
        end
        
        % 3)
        nsib1 = nmz;
        perm_ind = randsample(nsib,nsib1);
        jsib1 = jsib_all(perm_ind);
        isib1 = isib_all(perm_ind);
        
        %
        npair = nmz+nsib1+nur;
        
        jall = [jmz; jsib1; jur];
        iall = [imz; isib1; iur];
        
        all_label = cell(npair,1);
        all_label(1:nmz) = {'MZ'};
        all_label(nmz+1:nmz+nsib1) = {'SIB'};
        all_label(nmz+nsib1+1:npair) = {'UR'};
    end
    
    cv = cvpartition(all_label,'KFold',nfold);
    
    for ifold = 1:nfold
        train_ind = training(cv,ifold);
        test_ind = test(cv,ifold);
        
        [cvm_crit,acc,cmat,tpr,tnr,ppv,auc,sel_featnames] = run_fingerprint_train_test(dat,covar_dat,iall,jall,train_ind,test_ind,all_label,feat_names,analysis_flag,measure_flag,covar_flag,fs_flag,diff_func0,dist_func0,fitmodel,input_add);
        
        cvm_crit_all(itr,ifold) = cvm_crit;
        
        acc_all(itr,ifold) = acc;
        cmat_all(itr,ifold,:,:) = cmat;
        tpr_all(itr,ifold,:) = tpr;
        tnr_all(itr,ifold,:) = tnr;
        ppv_all(itr,ifold,:) = ppv;
        auc_all(itr,ifold,:) = auc;
        
        sel_feats_all(itr,ifold) = {sel_featnames};
    end
end

ML_STATS = [];
ML_STATS.cvm_acc_all = cvm_crit_all;
ML_STATS.acc_all = acc_all;
ML_STATS.cmat_all = cmat_all;
ML_STATS.tpr_all = tpr_all;
ML_STATS.tnr_all = tnr_all;
ML_STATS.ppv_all = ppv_all;
ML_STATS.auc_all = auc_all;
ML_STATS.sel_feats_all = sel_feats_all;

end
