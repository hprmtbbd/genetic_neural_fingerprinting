
clear,clc

workspace_path = 'D:\workspace';
proj_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript');
analysis_path = fullfile(proj_path,'h_analysis');

outdir_root = 'D:\Data\h_analysis_out';

func_path = fullfile(workspace_path,'MATLAB','functions');
func2_path = fullfile(proj_path,'mfile_h');

addpath(genpath(func_path),genpath(func2_path))

if isempty(gcp('nocreate')), parpool(16,'IdleTimeout',600), end

save_flag = 0;

data_flag = 1; % 1 - meg, 2 - fmri
measure_flag = 1; % 1 - global, 2 - local
analysis_flags = 1:4; % 1- three class, 2- MZ vs UR, 3- MZ vs SIB, 4- SIB vs UR

covar_flag = 0; % regress out covariates from features (graph measures)
fs_flag = 0;

diff_func = 'abs'; % abs corr

%% select MEG or fMRI
if measure_flag == 1
    measure_tag = 'GGM';
    Gnames = {'GE','CPL','T','S'};
elseif measure_flag == 2
    measure_tag = 'LGM';
    Gnames = {'STR','EVC','CC','NE'};
    nnode = 246;
end

if data_flag == 1
    infile = fullfile(analysis_path,['HCP_MEG_',measure_tag,'_Source_100%PropThres.csv']);
    data_tag = 'meg';
    
    connames = {'dwPLI','AEC','lcAEC'};
    freqBands = {'Delta','Theta','Alpha','lBeta','hBeta','lGamma'};
    nfreq = length(freqBands);
elseif data_flag == 2
    infile = fullfile(analysis_path,['HCP_fMRI_',measure_tag,'_Source_100%PropThres.csv']);
    data_tag = 'fmri';
    
    polarity_all = {'pos','neg'};
    connames = {'corr'};
    npol = length(polarity_all);
end

outfolder = 'ML_fingerprint_v2';
outfolder2 = ['featsel_',num2str(fs_flag)];
outdir = fullfile(outdir_root,outfolder,outfolder2,[measure_tag,'_',data_tag]);
if ~exist(outdir,'dir'),mkdir(outdir),end

%% load graph measures
nmethod = length(connames);
ngraph = length(Gnames);

%%%%
opts = detectImportOptions(infile);
in_data = readtable(infile,opts);
in_tit = in_data.Properties.VariableNames;
ggm_ind = 7:length(in_tit);

u_subjidlist = cellfun(@num2str,num2cell(in_data.('ID')),'UniformOutput',false);
ggm_dat1 = table2array(in_data(:,ggm_ind));
ggm_tit1 = in_tit(ggm_ind);
nsubj = length(u_subjidlist);

if measure_flag == 1
    ggm_dat = ggm_dat1;
    ggm_tit = ggm_tit1;
elseif measure_flag == 2
    sel_tit = contains(ggm_tit1,strcat('_',Gnames));
    ggm_dat = ggm_dat1(:,sel_tit);
    ggm_tit = ggm_tit1(sel_tit);
end

covar_ind = [3,4];
covar_dat = table2array(in_data(:,covar_ind));

if data_flag == 1
    m_ggm = []; % (isubj,ifr,igraph,imethod)
    m_tit = {};
    for ifr = 1:nfreq
        for igraph = 1:ngraph
            for imethod = 1:nmethod
                tit0 = [Gnames{igraph},'_',connames{imethod},'_',freqBands{ifr}];
                
                if measure_flag == 1
                    sel_ind = find(strcmp(tit0,ggm_tit));
                    ggm_dat0 = ggm_dat(:,sel_ind);
                    m_ggm(:,ifr,igraph,imethod) = ggm_dat0;
                elseif measure_flag == 2
                    sel_ind = find(contains(ggm_tit,tit0));
                    ggm_dat0 = ggm_dat(:,sel_ind);
                    m_ggm(:,:,ifr,igraph,imethod) = ggm_dat0;
                end
                
                m_tit(ifr,igraph,imethod) = {tit0};
            end
        end
    end
    
    nfeat = nfreq*ngraph;
elseif data_flag == 2
    nfeat = ngraph*npol;
    if measure_flag == 1
        m_ggm = ggm_dat; % (isubj,igraph*ipol)
        m_tit = ggm_tit;
    elseif measure_flag == 2
        m_ggm = reshape(ggm_dat,[nsubj,nnode,nfeat]);
        m_tit0 = reshape(ggm_tit,[nnode,nfeat]);
        m_tit = replace(m_tit0(1,:),'fMRI_Node1_','');
    end
end

%% construct relationship matrix
csv_path = analysis_path;
csv_infile = fullfile(csv_path,'HCP-YA_allSubjects_pedGT.csv');
opts = detectImportOptions(csv_infile);
opts.SelectedVariableNames = {'ID','FA','MO','FAMID','MZTWIN'};
csv_data = table2cell(readtable(csv_infile,opts));
csv_subjidlist = cellfun(@num2str,csv_data(:,1),'UniformOutput',false);
csv_faid = csv_data(:,2);
csv_moid = csv_data(:,3);
csv_famid = csv_data(:,4);
csv_mztwin = csv_data(:,5);

[~,loc] = ismember(u_subjidlist,csv_subjidlist);
faid = csv_faid(loc);
moid = csv_moid(loc);
famid = csv_famid(loc);
mztwin = csv_mztwin(loc);

phi = [];
for isubj = 1:nsubj
    subjid_i = u_subjidlist{isubj};
    faid_i = faid{isubj};
    moid_i = moid{isubj};
    mztwin_i = mztwin{isubj};
    
    for jsubj = 1:nsubj
        subjid_j = u_subjidlist{jsubj};
        faid_j = faid{jsubj};
        moid_j = moid{jsubj};
        mztwin_j = mztwin{jsubj};
        
        if isequal(subjid_i,subjid_j)
            phi0 = 1;
        elseif isequal(faid_i,faid_j) && isequal(moid_i,moid_j) && ~isempty(mztwin_i) && isequal(mztwin_i,mztwin_j)
            phi0 = 1;
        elseif isequal(faid_i,faid_j) && isequal(moid_i,moid_j)
            phi0 = 1/2;
        else
            phi0 = 0;
        end
        
        phi(isubj,jsubj) = phi0;
    end
end

phi0 = phi;
phi0(1:nsubj+1:end) = 0;

%%
phi0_liu = tril(phi0,-1);

[imz,jmz] = find(phi0_liu==1);

if data_flag == 1
   [isib,jsib] = find(phi0_liu==1/2);
   
   nmz = length(jmz); nsib = length(jsib);
   
   fam_locs00 = []; fam_locs = [];
elseif data_flag == 2
    u_famid = unique(famid,'stable');
    mztwin_lia = sum(phi0==1,2)~=0;
    
    fam_locs = {}; mz_lia_all = logical([]);
    for ifam = 1:length(u_famid)
        famid_i = u_famid{ifam};
        fam_ind = find(strcmp(famid_i,famid));
        mz_lia = any(mztwin_lia(fam_ind));
        fam_locs{ifam,1} = fam_ind';
        mz_lia_all(ifam,1) = mz_lia;
    end
    fam_size = cell2mat(cellfun(@length,fam_locs,'UniformOutput',0));
    
    fam_locs0 = fam_locs(~mz_lia_all);
    fam_size0 = fam_size(~mz_lia_all);
    fam_locs00 = fam_locs0(fam_size0 > 1);
    
    mz_fams = fam_locs(mz_lia_all);
    
    isib = cell2mat(cellfun(@(x) x(1),fam_locs00,'UniformOutput',0));
    jsib = cell2mat(cellfun(@(x) x(2),fam_locs00,'UniformOutput',0));
    
    nmz = length(jmz); nsib = length(jsib);
end

%%
niter = 1000; % 1000 100
nfold = 5;

col_tit = {'TPR','TNR','PPV','AUC'};
tit_a = {'MZ','SIB','UR'};
tits = {'MZ twins','Non-MZ siblings','Unrelated subjects','Weighted average'}';

for ia = 1:length(analysis_flags)
    analysis_flag = analysis_flags(ia);
    
    if ismember(analysis_flag,[2,3,4])
        two_class_flag = 1; % 1- classify only two classes
    else
        two_class_flag = 0;
    end
    
    %%
    if two_class_flag
        fitmodel = @fitcsvm;
    else
        fitmodel = @fitcecoc;
    end
    
    if isequal(fitmodel,@fitcsvm) || isequal(fitmodel,@fitcecoc)
        kfun = 'linear';
        C = 1;
        input_add = {'KernelFunction',kfun,'KernelScale','auto','BoxConstraint',C};
        
        if ~two_class_flag
            t = templateSVM(input_add{:});
            input_add = {'Learners',t};
        end
        model_tag = ['SVM-',kfun,'-',num2str(C)];
    end
    
    if covar_flag
        model_tag = [model_tag,'-','covar'];
    end
    
    if strcmp(diff_func,'abs')
        diff_func0 = @calc_abs_diff;
        dist_func0 = @calc_euclidean_dist;
    elseif strcmp(diff_func,'corr')
        diff_func0 = @calc_corr_diff;
        dist_func0 = @calc_corr_dist;
    end
    
    model_tag = [model_tag,'-',diff_func];
    
    %%
    if analysis_flag == 1
        row_tit = tits;
        class_tag = 'threeClass';
    elseif analysis_flag == 2
        row_tit = tits([1,3,4]);
        class_tag = [tit_a{1},'-',tit_a{3}];
    elseif analysis_flag == 3
        row_tit = tits([1,2,4]);
        class_tag = [tit_a{1},'-',tit_a{2}];
    elseif analysis_flag == 4
        row_tit = tits([2,3,4]);
        class_tag = [tit_a{2},'-',tit_a{3}];
    end
    
    outfile = fullfile(outdir,[model_tag,'_',class_tag,'.mat']);
    
    %%
    if ~exist(outfile,'file')
        
        ML_STATS_all = {};
        
        tic
        for imethod = 1:nmethod
            disp(connames{imethod})
            
            if data_flag == 1
                
                if measure_flag == 1
                    ggm0 = squeeze(m_ggm(:,:,:,imethod));
                    dat = reshape(ggm0,[nsubj,nfeat]);
                    feat_names = reshape(m_tit(:,:,imethod),[nfeat,1]);
                elseif measure_flag == 2
                    ggm0 = squeeze(m_ggm(:,:,:,:,imethod));
                    dat = reshape(ggm0,nsubj,[],nfeat);
                    feat_names = reshape(m_tit(:,:,imethod),[nfeat,1]);
                end
                
            elseif data_flag == 2
                dat = m_ggm;
                feat_names = m_tit;
            end
            
            ML_STATS = run_fingerprint(dat,covar_dat,feat_names,imz,jmz,isib,jsib,phi,fam_locs,fam_locs00,nsubj,nmz,nsib,analysis_flag,data_flag,measure_flag,covar_flag,nfold,niter,fs_flag,diff_func0,dist_func0,fitmodel,input_add);
            
            ML_STATS_all{imethod,1} = ML_STATS;
        end
        toc
        
        %%
        if save_flag
            save(outfile,'ML_STATS_all')
        end
        
    else
        load(outfile)
    end
    
    %%
    disp('-----------------------------------------------------------------')
    disp(class_tag)
    c_all = {};
    for imethod = 1:nmethod
        ML_STATS = ML_STATS_all{imethod};
        
        conn_tit = connames{imethod};
        con_tits = cell(1,4);
        con_tits(1) = compose('%s (acc = %0.3f +/ %0.3f)',conn_tit,squeeze(nanmean(ML_STATS.acc_all,[1,2])),squeeze(nanstd(ML_STATS.acc_all,[],[1,2])));
        disp(' ')
        disp(conn_tit)
        disp(['ACC: ',num2str(squeeze(nanmean(ML_STATS.acc_all,[1,2])))])
        disp(['TPR: ',num2str(squeeze(nanmean(ML_STATS.tpr_all,[1,2]))')])
        disp(['TNR: ',num2str(squeeze(nanmean(ML_STATS.tnr_all,[1,2]))')])
        disp(['PPV: ',num2str(squeeze(nanmean(ML_STATS.ppv_all,[1,2]))')])
        disp(['AUC: ',num2str(squeeze(nanmean(ML_STATS.auc_all,[1,2]))')])
        
        c1 = compose('%0.3f +/- %0.3f',squeeze(nanmean(ML_STATS.tpr_all,[1,2])),squeeze(nanstd(ML_STATS.tpr_all,[],[1,2])));
        c2 = compose('%0.3f +/- %0.3f',squeeze(nanmean(ML_STATS.tnr_all,[1,2])),squeeze(nanstd(ML_STATS.tnr_all,[],[1,2])));
        c3 = compose('%0.3f +/- %0.3f',squeeze(nanmean(ML_STATS.ppv_all,[1,2])),squeeze(nanstd(ML_STATS.ppv_all,[],[1,2])));
        c4 = compose('%0.3f +/- %0.3f',squeeze(nanmean(ML_STATS.auc_all,[1,2])),squeeze(nanstd(ML_STATS.auc_all,[],[1,2])));
        cc = [con_tits; col_tit; c1,c2,c3,c4];
        
        c_all = [c_all; [[{[]; []}; row_tit],cc]];
    end
    disp('-----------------------------------------------------------------')
    
    if save_flag
        outfile_xls = fullfile(outdir,[model_tag,'_',class_tag,'.xlsx']);
        writecell(c_all,outfile_xls)
    end
    
end

%%

function d_scores = calc_abs_diff(scores1,scores2)

d_scores = abs(scores1-scores2);

end

function d_scores = calc_corr_diff(scores1,scores2)

c_scores1 = scores1-mean(scores1,2);
c_scores2 = scores2-mean(scores2,2);
d_scores = 1 - c_scores1.*c_scores2./sqrt(sum(c_scores1.^2,2).*sum(c_scores2.^2,2));

end

%%
function d_scores = calc_euclidean_dist(scores1,scores2)

d_scores = squeeze(sqrt(sum((scores1-scores2).^2,2)));

end

function d_scores = calc_corr_dist(scores1,scores2)

c_scores1 = scores1-mean(scores1,2);
c_scores2 = scores2-mean(scores2,2);
d_scores = squeeze(1 - sum(c_scores1.*c_scores2,2)./sqrt(sum(c_scores1.^2,2).*sum(c_scores2.^2,2)));

end
