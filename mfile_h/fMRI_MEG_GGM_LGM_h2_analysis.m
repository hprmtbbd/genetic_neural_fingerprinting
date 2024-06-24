
clear,clc

indir_root = 'E:'; % F
indir = fullfile(indir_root,'UT-Austin','h_analysis');
indir_meg = fullfile(indir,'phaatef_MEG_LGM_genetictest');
indir_fmri = fullfile(indir,'phaatef_fMRI_LGM_genetictest');
indir_meg_ggm = fullfile(indir,'phaatef_MEG_GGM_genetictest');
indir_fmri_ggm = fullfile(indir,'phaatef_fMRI_GGM_genetictest_v2');

workspace_path = 'D:\workspace';
toolbox_path = fullfile(workspace_path,'MATLAB','toolbox');
bnet_path = fullfile(toolbox_path,'BrainNetViewer_20191031');
func_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript','mfile_h');
func_path2 = fullfile(workspace_path,'MATLAB','functions');
addpath(genpath(bnet_path),genpath(func_path),genpath(func_path2))

proj_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript');
analysis_path = fullfile(proj_path,'h_analysis');
matfile_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript','mfile_h','HCP');

outpath_root = 'D:\Data\h_analysis_out';

outdir_fig = fullfile(proj_path,'Figures','version3');

nodelabel_file = fullfile(matfile_path,'brainnetomeLabel.mat');
nodeloc_file = fullfile(matfile_path,'brainnetomeMedoidLocation_4mm.mat');
load(nodelabel_file), load(nodeloc_file)
nnode = length(nodelabel);

screenPos = get(0,'Screensize');

save_flag_png = 0;
save_flag_fig = 0;

bnet_save_flag = 0;
xlsout_save_flag = 0;

if ~exist(outdir_fig,'dir'),mkdir(outdir_fig),end

%% load meg LGM h2 values
outfile = fullfile(outpath_root,'SOLAR','MEG_LGM_Source_100%PropThres_hval.mat');

if ~exist(outfile,'file')
    infile_LGM = fullfile(indir_meg,'HCP_MEG_LGM_Source_100%PropThres.csv');
    tab = readtable(infile_LGM);
    
    tit = tab.Properties.VariableNames;
    tit_lgm = tit(7:end);
    
    ind_CoC = contains(tit_lgm,'_CoC');
    tit_lgm(ind_CoC) = [];
    
    conmethods = {'dwPLI','AEC','lcAEC'};
    freqBands = {'Delta','Theta','Alpha','lBeta','hBeta','lGamma'};
    graphMeasures = {'STR','EVC','CC','NE'};
    
    h_all = []; p_all = []; he_all = []; cov_all = []; kurt_all = []; kurt_score_all = {};
    parfor ilgm = 1:length(tit_lgm)
        varname = tit_lgm{ilgm};
        splt_varname = strsplit(varname,'_');
        gm_name = splt_varname{2};
        
        indir0 = fullfile(indir_meg,gm_name,varname);
        [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0);
        
        h_all(ilgm,1) = hval;
        p_all(ilgm,1) = pval;
        he_all(ilgm,1) = herror;
        cov_all(ilgm,1) = covar_var;
        kurt_all(ilgm,1) = kurt;
        kurt_score_all{ilgm,1} = kurt_score;
    end
    
    pFDR_all = mafdr(p_all,'BHFDR',1);
    
    MEG_LGM_H2 = [];
    MEG_LGM_H2.h_all = h_all;
    MEG_LGM_H2.p_all = p_all;
    MEG_LGM_H2.pFDR_all = pFDR_all;
    MEG_LGM_H2.he_all = he_all;
    MEG_LGM_H2.cov_all = cov_all;
    MEG_LGM_H2.kurt_all = kurt_all;
    MEG_LGM_H2.kurt_score_all = kurt_score_all;
    MEG_LGM_H2.tit_lgm = tit_lgm;
    MEG_LGM_H2.conmethods = conmethods;
    MEG_LGM_H2.freqBands = freqBands;
    MEG_LGM_H2.graphMeasures = graphMeasures;
    
    save(outfile,'MEG_LGM_H2')
else
    load(outfile)
end

%% load fmri LGM h2 values
outfile = fullfile(outpath_root,'SOLAR_fMRI','fMRI_LGM_100%PropThres_hval.mat');

if ~exist(outfile,'file')
    infile_LGM = fullfile(analysis_path,'HCP_fMRI_LGM_Source_100%PropThres.csv');
    T = readtable(infile_LGM);
    tit = T.Properties.VariableNames;
    fmri_subjid = cellfun(@num2str,num2cell(T.ID),'UniformOutput',false);
    
    tit_lgm = tit(7:end);
    
    ind_CoC = contains(tit_lgm,'_CoC');
    tit_lgm(ind_CoC) = [];
    
    graphMeasures = {'STR','EVC','CC','NE'};
    
    polarity_all = {'pos','neg'};
    conmethods = {'corr'};
    
    h_all = []; p_all = []; he_all = []; cov_all = []; kurt_all = []; kurt_score_all = {};
    for iggm = 1:length(tit_lgm)
        varname0 = tit_lgm{iggm};
        tok = strsplit(varname0,'_');
        gm_name = tok{3};
        varname = erase(varname0,'fMRI_');
        
        indir0 = fullfile(indir_fmri,gm_name,varname);
        [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0);
        
        h_all(iggm,1) = hval;
        p_all(iggm,1) = pval;
        he_all(iggm,1) = herror;
        cov_all(iggm,1) = covar_var;
        kurt_all(iggm,1) = kurt;
        kurt_score_all{iggm,1} = kurt_score;
    end
    
    pFDR_all = mafdr(p_all,'BHFDR',1);
    
    fMRI_LGM_H2 = [];
    fMRI_LGM_H2.h_all = h_all;
    fMRI_LGM_H2.p_all = p_all;
    fMRI_LGM_H2.pFDR_all = pFDR_all;
    fMRI_LGM_H2.he_all = he_all;
    fMRI_LGM_H2.cov_all = cov_all;
    fMRI_LGM_H2.kurt_all = kurt_all;
    fMRI_LGM_H2.kurt_score_all = kurt_score_all;
    fMRI_LGM_H2.tit_lgm = tit_lgm;
    fMRI_LGM_H2.conmethods = conmethods;
    fMRI_LGM_H2.polarity_all = polarity_all;
    fMRI_LGM_H2.graphMeasures = graphMeasures;
    
    save(outfile,'fMRI_LGM_H2')
else
    load(outfile)
end

%% load meg GGM h2 values
outfile = fullfile(outpath_root,'SOLAR','MEG_GGM_Source_100%PropThres_hval.mat');

if ~exist(outfile,'file')
    infile_GGM = fullfile(indir_meg_ggm,'HCP_MEG_GGM_Source_100%PropThres.csv');
    tab = readtable(infile_GGM);
    
    tit = tab.Properties.VariableNames;
    tit_ggm = tit(7:end);
    
    conmethods = {'dwPLI','AEC','lcAEC'};
    freqBands = {'Delta','Theta','Alpha','lBeta','hBeta','lGamma'};
    graphMeasures = {'GE','CPL','T','S'};
    
    h_all = []; p_all = []; he_all = []; cov_all = []; kurt_all = []; kurt_score_all = {};
    for iggm = 1:length(tit_ggm)
        varname = tit_ggm{iggm};
        indir0 = fullfile(indir_meg_ggm,varname);
        
        [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0);
        
        h_all(iggm,1) = hval;
        p_all(iggm,1) = pval;
        he_all(iggm,1) = herror;
        cov_all(iggm,1) = covar_var;
        kurt_all(iggm,1) = kurt;
        kurt_score_all{iggm,1} = kurt_score;
    end
    
    pFDR_all = mafdr(p_all,'BHFDR',1);
    
    MEG_GGM_H2 = [];
    MEG_GGM_H2.h_all = h_all;
    MEG_GGM_H2.p_all = p_all;
    MEG_GGM_H2.pFDR_all = pFDR_all;
    MEG_GGM_H2.he_all = he_all;
    MEG_GGM_H2.cov_all = cov_all;
    MEG_GGM_H2.kurt_all = kurt_all;
    MEG_GGM_H2.kurt_score_all = kurt_score_all;
    MEG_GGM_H2.tit_ggm = tit_ggm;
    MEG_GGM_H2.conmethods = conmethods;
    MEG_GGM_H2.freqBands = freqBands;
    MEG_GGM_H2.graphMeasures = graphMeasures;
    
    save(outfile,'MEG_GGM_H2')
else
    load(outfile)
end

%% load fmri GGM h2 values
outfile = fullfile(outpath_root,'SOLAR_fMRI','fMRI_GGM_100%PropThres_hval.mat');

if ~exist(outfile,'file')
    infile_GGM = fullfile(indir_fmri_ggm,'HCP_fMRI_GGM_Source_100%PropThres.csv');
    tab = readtable(infile_GGM);
    
    tit = tab.Properties.VariableNames;
    tit_ggm = tit(7:end);
    
    polarity_all = {'pos','neg'};
    conmethods = {'corr'};
    graphMeasures = {'GE','CPL','T','S'};
    
    h_all = []; p_all = []; he_all = []; cov_all = []; kurt_all = []; kurt_score_all = {};
    for iggm = 1:length(tit_ggm)
        varname0 = tit_ggm{iggm};
        varname = erase(varname0,'fMRI_');
        
        indir0 = fullfile(indir_fmri_ggm,varname);
        [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0);
        
        h_all(iggm,1) = hval;
        p_all(iggm,1) = pval;
        he_all(iggm,1) = herror;
        cov_all(iggm,1) = covar_var;
        kurt_all(iggm,1) = kurt;
        kurt_score_all{iggm,1} = kurt_score;
    end
    
    pFDR_all = mafdr(p_all,'BHFDR',1);
    
    fMRI_GGM_H2 = [];
    fMRI_GGM_H2.h_all = h_all;
    fMRI_GGM_H2.p_all = p_all;
    fMRI_GGM_H2.pFDR_all = pFDR_all;
    fMRI_GGM_H2.he_all = he_all;
    fMRI_GGM_H2.cov_all = cov_all;
    fMRI_GGM_H2.kurt_all = kurt_all;
    fMRI_GGM_H2.kurt_score_all = kurt_score_all;
    fMRI_GGM_H2.tit_ggm = tit_ggm;
    fMRI_GGM_H2.conmethods = conmethods;
    fMRI_GGM_H2.polarity_all = polarity_all;
    fMRI_GGM_H2.graphMeasures = graphMeasures;
    
    save(outfile,'fMRI_GGM_H2')
else
    load(outfile)
end

%% compute FDR-corrected p-values for GGM h2
p_GGM = [MEG_GGM_H2.p_all; fMRI_GGM_H2.p_all];
pFDR_GGM = mafdr(p_GGM,'BH',1);
tit_ggm = [strcat('MEG_',MEG_GGM_H2.tit_ggm'); fMRI_GGM_H2.tit_ggm'];
sig_ind = pFDR_GGM<0.05;
sig_tit = tit_ggm(sig_ind);

pFDR_GGM_MEG = pFDR_GGM(contains(tit_ggm,'MEG_'));
pFDR_GGM_fMRI = pFDR_GGM(contains(tit_ggm,'fMRI_'));

GGM_names = {'Global Efficiency','Characteristic Path Length','Transitivity','Synchronizability'};
polarity_names = {'Positive Correlations','Negative Correlations'};

%% make MEG and fMRI GGM h2 figure

figPos_a = [50,100,screenPos(3)*0.9,screenPos(4)*0.5];

nmethod = length(MEG_GGM_H2.conmethods);
nfreq = length(MEG_GGM_H2.freqBands);

npol = length(fMRI_GGM_H2.polarity_all);
ngraph = length(fMRI_GGM_H2.graphMeasures);

pFDR_mat1 = reshape(pFDR_GGM_fMRI,[ngraph,npol]);
h_mat1 = reshape(fMRI_GGM_H2.h_all,[ngraph,npol]);
he_mat1 = reshape(fMRI_GGM_H2.he_all,[ngraph,npol]);

pFDR_mat2 = reshape(pFDR_GGM_MEG,[ngraph,nmethod,nfreq]);
h_mat2 = reshape(MEG_GGM_H2.h_all,[ngraph,nmethod,nfreq]);
he_mat2 = reshape(MEG_GGM_H2.he_all,[ngraph,nmethod,nfreq]);

FDR_flag = 1;
alpha = [0.05,0.005];

for igraph = 1:ngraph
    GE_name = fMRI_GGM_H2.graphMeasures{igraph};
    
    %-------------------------------------------------------
    figure_name = ['MEG_fMRI_',GE_name,'_hval'];
    if FDR_flag
        figure_name = [figure_name,'_FDR'];
    end
    
    h = figure;
    set(h,'name',figure_name,'position',figPos_a)
    
    %%%%
    % (ifr,imethod)
    p_values = squeeze(pFDR_mat2(igraph,:,:))';
    m = squeeze(h_mat2(igraph,:,:))';
    se = squeeze(he_mat2(igraph,:,:))'; se(isnan(se)) = 0;
    
    subplot(1,2,1)
    mu = m;
    sigma = se;
    p_val = p_values;
    
    C = [0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880];
    CC = permute(C,[3,1,2]);
    
    ax = superbar(mu,'BarFaceColor',CC,'E',sigma,'ErrorbarColor','k','P',p_val,...
        'PStarShowGT',false,'PStarThreshold',alpha,'PStarIcon','*',...
        'PStarShowNS',false,'PStarFontSize',14,'PLineColor','k','PLineWidth',2,'PStarColor','k',...
        'PStarOffset',0.02);
    
    [m1,n1] = size(mu');
    loc1 = 0.15* linspace(-ceil(m1/2),ceil(m1/2),m1)'*ones(1,n1) + ones(m1,1)*(1:n1);
    xlim( [ 0.5 ( max(loc1(:))+min(loc1(:))/2 ) ] )
    xticks(1:nfreq)
    set(gca,'XTickLabel',MEG_GGM_H2.freqBands,'FontName','Arial')
    set(gca,'Fontsize',16)
    set(gca,'YGrid','on')
    ylim([-0.15,1]),yticks(0:0.1:1)
    ylabel('Heritability (h^2)')
    tsub = title('MEG');
    set(tsub,'FontSize',18)
    
    fs = gca;
    pos_s = get(fs,'Position');
    newpos_s = [pos_s(1) pos_s(2)-0.025 pos_s(3)*1.3 pos_s(4)*0.95];
    set(fs,'Position',newpos_s)
    
    lgd = legend(ax(1,:),MEG_GGM_H2.conmethods,'Orientation','Vertical');
    
    %%%%
    p_values = pFDR_mat1(igraph,:);
    m = h_mat1(igraph,:);
    se = he_mat1(igraph,:);
    
    subplot(1,2,2)
    mu = repmat(m,[2,1]);
    sigma = repmat(se,[2,1]);
    p_val = repmat(p_values,[2,1]);
    
    C = [1 0 0; 0.3010 0.7450 0.9330];
    
    [HB, HE, HPT] = superbar(mu,'BarFaceColor',permute(C,[3,1,2]),'E',sigma,'ErrorbarColor','k','P',p_val,...
        'PStarShowGT',false,'PStarThreshold',alpha,'PStarIcon','*',...
        'PStarShowNS',false,'PStarFontSize',14,'PLineColor','k','PLineWidth',2,'PStarColor','k',...
        'PStarOffset',0.02);
    
    delete(HB(2,1)), delete(HB(2,2))
    delete(HE(2,1)), delete(HE(2,2))
    delete(HPT(2,1)), delete(HPT(2,2))
    
    xticks(1)
    set(gca,'XTickLabel','Pearson Correlation')
    set(gca,'Fontsize',16)
    set(gca,'YGrid','on')
    ylim([-0.15,1]),yticks(0:0.1:1)
    ylabel('Heritability (h^2)')
    tsub = title('fMRI');
    set(tsub,'FontSize',18)
    
    fs = gca;
    pos_s = get(fs,'Position');
    newpos_s = [pos_s(1)+0.075 pos_s(2)-0.025 pos_s(3)*0.55 pos_s(4)*0.95];
    set(fs,'Position',newpos_s)
    
    lgd = legend(polarity_names,'Orientation','Vertical');
    pos_l = get(lgd,'Position');
    
    sgt = annotation('textbox','String',[GGM_names{igraph},' (',GE_name,')'],'EdgeColor','none','FontSize',18);
    set(sgt,'FitBoxToText','on')
    drawnow
    tl_pos = get(sgt,'Position'); xtl = 0.5-tl_pos(3)/2+0.05;
    set(sgt,'Position',[xtl 1 1 0])
    
    dpi = '1200';
    img_outfile = [outdir_fig,'\',figure_name];
    img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(h,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
end

%% meg vs fmri LGM h2 (Dice coefficient, signed-rank test)
nfeat_meg = length(MEG_LGM_H2.graphMeasures)*length(MEG_LGM_H2.freqBands)*length(MEG_LGM_H2.conmethods);
meg_tit_lgm = reshape(MEG_LGM_H2.tit_lgm',[nnode,nfeat_meg]);
meg_feat = replace(meg_tit_lgm(1,:),'Node1_','MEG_')';
h_mat_meg = reshape(MEG_LGM_H2.h_all,[nnode,nfeat_meg]);

LGM_names = {'Strength','Eigenvector Centrality','Clustering Coefficient','Nodal Efficiency'};

%%%%
nfeat_fmri = length(fMRI_LGM_H2.graphMeasures)*length(fMRI_LGM_H2.polarity_all)*length(fMRI_LGM_H2.conmethods);
fmri_tit_lgm = reshape(fMRI_LGM_H2.tit_lgm',[nnode,nfeat_fmri]);
fmri_feat = replace(fmri_tit_lgm(1,:),'Node1_','')';
h_mat_fmri = reshape(fMRI_LGM_H2.h_all,[nnode,nfeat_fmri]);

feat = [meg_feat; fmri_feat];
h_mat = [h_mat_meg,h_mat_fmri];

%%%% compute FDR-corrected p-values for LGM h2
p_LGM = [MEG_LGM_H2.p_all; fMRI_LGM_H2.p_all];
pFDR_LGM = mafdr(p_LGM,'BH',1);
tit_lgm = [strcat('MEG_',MEG_LGM_H2.tit_lgm'); fMRI_LGM_H2.tit_lgm'];
sig_lgm_ind = pFDR_LGM<0.05;
sig_lgm_tit = tit_lgm(sig_lgm_ind);

pFDR_LGM_MEG = pFDR_LGM(contains(tit_lgm,'MEG_'));
pFDR_LGM_fMRI = pFDR_LGM(contains(tit_lgm,'fMRI_'));

pFDR_mat_meg = reshape(pFDR_LGM_MEG,[nnode,nfeat_meg]);
pFDR_mat_fmri = reshape(pFDR_LGM_fMRI,[nnode,nfeat_fmri]);

num_sigFDR_meg = sum(pFDR_mat_meg < 0.05,1)'/nnode;
num_sigFDR_fmri = sum(pFDR_mat_fmri < 0.05,1)'/nnode;

%%%%
DC_thres = 0.5;

DC_fmri = []; TT_fmri = []; TP_fmri = []; DT_fmri = [];
for ifeat = 1:length(fmri_feat)
    fmri_feat0 = fmri_feat{ifeat};
    tok = strsplit(fmri_feat0,'_');
    gm_name = tok{2};
    meg_gm_ind = contains(meg_feat,gm_name);
    meg_gm_feat = meg_feat(meg_gm_ind);
    
    h_mat_fmri0 = h_mat_fmri(:,ifeat);
    h_mat_meg0 = h_mat_meg(:,meg_gm_ind);
    
    pFDR_mat_fmri0 = pFDR_mat_fmri(:,ifeat);
    sig_mat_fmri0 = pFDR_mat_fmri0 < 0.05;
    high_h_fmri0 = h_mat_fmri0 >= quantile(h_mat_fmri0,DC_thres);
    lia_mat_fmri0 = sig_mat_fmri0 & high_h_fmri0;
    
    pFDR_mat_meg0 = pFDR_mat_meg(:,meg_gm_ind);
    sig_mat_meg0 = pFDR_mat_meg0 < 0.05;
    high_h_meg0 = h_mat_meg0 >= quantile(h_mat_meg0,DC_thres);
    lia_mat_meg0 = sig_mat_meg0 & high_h_meg0;
    
    nC = sum((lia_mat_fmri0 & lia_mat_meg0),1);
    nT = (sum(lia_mat_fmri0) + sum(lia_mat_meg0,1));
    dc = 2*nC./nT;
    
    ptt = []; tt = []; dt = [];
    for ifeat_meg = 1:length(meg_gm_feat)
        [pt,ht,stat_t] = signrank(h_mat_fmri0,h_mat_meg0(:,ifeat_meg));
        ptt(ifeat_meg,1) = pt;
        tt(ifeat_meg,1) = stat_t.zval;
        dt(ifeat_meg,1) =  mean(h_mat_meg0(:,ifeat_meg)-h_mat_fmri0);
    end
    
    DC_fmri(ifeat,:) = dc;
    
    TT_fmri(ifeat,:) = tt;
    PT_fmri(ifeat,:) = ptt;
    DT_fmri(ifeat,:) = dt;
end

PT_fmri_FDR_rs = mafdr(PT_fmri(:),'BH',1);
PT_fmri_FDR = reshape(PT_fmri_FDR_rs,size(PT_fmri));
PT_sig = PT_fmri_FDR < 0.05;

tok1 = cellfun(@(x) strsplit(x,'_'),meg_feat,'UniformOutput',0);
u_meg_feat = replace(unique(cellfun(@(x) [x{end-1},'_',x{end}],tok1,'UniformOutput',0),'stable'),'_',' ');

u_fmri_feat = replace(fmri_feat,'fMRI_','');
fmri_feat_tit = replace(u_fmri_feat,{'corr_pos','corr_neg','_'},{'Positive Correlations','Negative Correlations',' '});
reorder_fmri = [1,ngraph+1,2,ngraph+2,3,ngraph+3,4,ngraph+4];

UT_fmri = DT_fmri;
UT_fmri(~PT_sig) = NaN;
for imethod = 1:length(MEG_LGM_H2.conmethods)
    conmethod = MEG_LGM_H2.conmethods{imethod};
    sel_ind = startsWith(u_meg_feat,[conmethod,' ']);
    
    ytL = fmri_feat_tit(reorder_fmri);
    xtL = u_meg_feat(sel_ind);
    yytL = fMRI_LGM_H2.graphMeasures;
    xlab = 'MEG Local Graph Measures';
    ylab = 'fMRI Local Graph Measures';
    yylab = '';
    
    % dice coefficient = sig heritable brain regions have high overlap between MEG and fMRI
    dat = DC_fmri(reorder_fmri,sel_ind);
    clims = [0,round(max(DC_fmri(:))*10)/10];
    cmap = brewermap([],'Reds');
    tit_lab = ['Dice similarity coefficient (DSC) between fMRI and MEG ',conmethod];
    
    fig = figure('Position',[50,100,screenPos(3)*0.5,screenPos(4)*0.6]);
    yyplot_matrix(dat,xtL,ytL,yytL,xlab,ylab,yylab,tit_lab,clims,cmap);
    
    figure_name = ['MEG_fMRI_DSC_',conmethod,'_hval'];
    dpi = '1200';
    img_outfile = [outdir_fig,'\',figure_name];
    img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
    
    % sign rank test - MEG vs fMRI
    dat = UT_fmri(reorder_fmri,sel_ind); 
    clims = [-1,1]*max(abs(dat(:)));
    cmap = flip(brewermap([],'RdBu'));
    tit_lab = ['Signed-rank test for h^2(MEG ',conmethod,')-h^2(fMRI)'];
    
    fig = figure('Position',[50,100,screenPos(3)*0.5,screenPos(4)*0.6]);
    yyplot_matrix(dat,xtL,ytL,yytL,xlab,ylab,yylab,tit_lab,clims,cmap);
    
    figure_name = ['MEG_fMRI_ttest_',conmethod,'_hval'];
    dpi = '1200';
    img_outfile = [outdir_fig,'\',figure_name];
    img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
end

% colorbars
figure_name = ['colorbar1'];
fig = figure;
imagesc(checkerboard(10),[0,round(max(DC_fmri(:))*10)/10])
colormap(brewermap([],'Reds'))
colorbar
dpi = '1200';
img_outfile = [outdir_fig,'\',figure_name];
img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)

figure_name = ['colorbar2'];
fig = figure;
imagesc(checkerboard(10),[-1,1]*max(abs(DT_fmri(:))))
colormap(flip(brewermap([],'RdBu')))
colorbar
dpi = '1200';
img_outfile = [outdir_fig,'\',figure_name];
img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)

%% plot MEG LGM h2 (BrainNet Viewer)
outpath = fullfile(outpath_root,'SOLAR');
cfg_file = fullfile(outpath,'bnet_config2.mat');
surface_file = fullfile(bnet_path,'Data','SurfTemplate','BrainMesh_ICBM152.nv');

tab_out = {};
tab2_out = {};
for imethod = 1:length(MEG_LGM_H2.conmethods)
    conmethod = MEG_LGM_H2.conmethods{imethod};
    ind_method = contains(MEG_LGM_H2.tit_lgm,['_',conmethod,'_']);
    
    con_tit = {}; con_tit(1:length(MEG_LGM_H2.freqBands)+1) = {''}; con_tit(1) = {conmethod};
    
    vec_sq = MEG_LGM_H2.graphMeasures';
    vec2_sq = MEG_LGM_H2.graphMeasures';
    for ifr = 1:length(MEG_LGM_H2.freqBands)
        freqBand = MEG_LGM_H2.freqBands{ifr};
        ind_freq = contains(MEG_LGM_H2.tit_lgm,freqBand);
        
        vec = {};
        vec2 = {};
        for igraph = 1:length(MEG_LGM_H2.graphMeasures)
            Gname = MEG_LGM_H2.graphMeasures{igraph};
            ind_gname = contains(MEG_LGM_H2.tit_lgm,['_',Gname,'_']);
            
            tit0 = [Gname,'_',conmethod,'_',freqBand];
            
            png_outdir = fullfile(outpath,Gname,'bnet_png_FDR');
            if ~exist(png_outdir,'dir')
                mkdir(png_outdir)
            end
            outdir_node_BNV = fullfile(outpath,Gname,'node_files');
            if ~exist(outdir_node_BNV,'dir')
                mkdir(outdir_node_BNV)
            end
            
            bnet_outtit = ['HCP_MEG_100%PropThres_',tit0,'_hvalGT_FDR'];
            img_outfile = fullfile(png_outdir,[bnet_outtit,'.png']);
            
            ind = ind_gname & ind_method & ind_freq;
            tit_lgm0 = MEG_LGM_H2.tit_lgm(ind)';
            h = MEG_LGM_H2.h_all(ind);
            he = MEG_LGM_H2.he_all(ind);
            pFDR = pFDR_LGM_MEG(ind);
            zt = abs(norminv(pFDR./2));
            
            hmean = mean(h);
            hstd = std(h);
            
            hmed = median(h);
            h1 = quantile(h,0.25);
            h2 = quantile(h,0.75);
            % entry2 = compose('%0.3f (%0.3f-%0.3f)',hmed,h1,h2);
            entry2 = compose('%0.3f (%0.3f)',hmean,hstd);
            vec2 = [vec2; entry2];
            
            num_sig = length(find(pFDR<0.05));
            num_zero = length(find(h==0));
            entry = compose('[%0.3f (%0.3f), %d, %d]',hmean,hstd,num_sig,num_zero);
            vec = [vec; entry];
            
            blank = {}; blank(1:length(nodelabel),1) = {'-'};
            node_data = [num2cell(nodeloc),num2cell(h),num2cell(zt),blank];
            outfile_node = fullfile(outdir_node_BNV,[bnet_outtit,'.node']);
            if ~exist(outfile_node,'file')
                node_fid = fopen(outfile_node,'a+');
                for i = 1:length(node_data)
                    fprintf(node_fid,'%.3f %.3f %.3f %f %f %s\n',node_data{i,:});
                end
                fclose(node_fid);
            end
            
            bnet_input = {surface_file,outfile_node,cfg_file};
            if bnet_save_flag
                bnet_input{end+1} = img_outfile;
                BrainNet_MapCfg(bnet_input{:});
            end
        end
        
        vec_sq = [vec_sq,vec];
        vec2_sq = [vec2_sq,vec2];
    end
    
    vec_sqq = [[{''},MEG_LGM_H2.freqBands]; vec_sq];
    vec2_sqq = [[{''},MEG_LGM_H2.freqBands]; vec2_sq];
    tab_out = [tab_out; con_tit; vec_sqq];
    tab2_out = [tab2_out; con_tit; vec2_sqq];
end

outfile = fullfile(outdir_fig,'HCP_MEG_LGM_100%PropThres_hvalGT.xlsx');
if xlsout_save_flag
    xlswrite(outfile,tab_out,'SOLAR')
end

%% plot fMRI LGM h2 (BrainNet Viewer)
outpath = fullfile(outpath_root,'SOLAR_fMRI');
cfg_file = fullfile(outpath,'bnet_config2.mat');
surface_file = fullfile(bnet_path,'Data','SurfTemplate','BrainMesh_ICBM152.nv');

conname = fMRI_LGM_H2.conmethods{1};

vec_sq = fMRI_LGM_H2.graphMeasures';
for ipol = 1:length(fMRI_LGM_H2.polarity_all)
    polarity = fMRI_LGM_H2.polarity_all{ipol};
    
    vec = {};
    for igraph = 1:length(fMRI_LGM_H2.graphMeasures)
        Gname = fMRI_LGM_H2.graphMeasures{igraph};
        
        tit0 = [Gname,'_',conname,'_',polarity];
        tit_ind = contains(fMRI_LGM_H2.tit_lgm,tit0);
        tit_gm0 = fMRI_LGM_H2.tit_lgm(tit_ind)';
        
        png_outdir = fullfile(outpath,Gname,'bnet_png_FDR');
        if ~exist(png_outdir,'dir')
            mkdir(png_outdir)
        end
        outdir_node_BNV = fullfile(outpath,Gname,'node_files');
        if ~exist(outdir_node_BNV,'dir')
            mkdir(outdir_node_BNV)
        end
        
        bnet_outtit = ['HCP_fMRI_100%PropThres_',tit0,'_hvalGT_FDR'];
        img_outfile = fullfile(png_outdir,[bnet_outtit,'.png']);
        
        h = fMRI_LGM_H2.h_all(tit_ind);
        p_FDR = pFDR_LGM_fMRI(tit_ind);
        zt = abs(norminv(p_FDR./2));
        
        hmean = mean(h);
        hstd = std(h);
        num_sig = length(find(p_FDR<0.05));
        num_zero = length(find(h==0));
        entry = compose('[%0.3f (%0.3f), %d, %d]',hmean,hstd,num_sig,num_zero);
        vec = [vec; entry];
        
        blank = {}; blank(1:length(nodelabel),1) = {'-'};
        node_data = [num2cell(nodeloc),num2cell(h),num2cell(zt),blank];
        outfile_node = fullfile(outdir_node_BNV,[bnet_outtit,'.node']);
        if ~exist(outfile_node,'file')
            node_fid = fopen(outfile_node,'a+');
            for i = 1:length(node_data)
                fprintf(node_fid,'%.3f %.3f %.3f %f %f %s\n',node_data{i,:});
            end
            fclose(node_fid);
        end
        
        bnet_input = {surface_file,outfile_node,cfg_file};
        if bnet_save_flag
            bnet_input{end+1} = img_outfile;
            BrainNet_MapCfg(bnet_input{:});
        end
    end
    vec_sq = [vec_sq,vec];
    
end
tab_out = [[{''},fMRI_LGM_H2.polarity_all]; vec_sq];

outfile = fullfile(outdir_fig,'HCP_fMRI_LGM_100%PropThres_hvalGT.xlsx');
if xlsout_save_flag
    xlswrite(outfile,tab_out,'SOLAR')
end

%% k-means clustering of LGM h2 and make tables of summary stats for each cluster
distance_metric = 'sqeuclidean';
options = statset('UseParallel',1);
clust_func = @(X,k) kmeans(X,k,'Distance',distance_metric,'MaxIter',1e3,'Replicates',1e3,'Options',options);

mh_mat_meg = mean(h_mat_meg,1)';
sh_mat_meg = std(h_mat_meg,[],1)';

p_mat_meg = reshape(MEG_LGM_H2.pFDR_all,[nnode,nfeat_meg]);
group_p_mat = mean(p_mat_meg,1)';
sig_feat_ind = true(size(meg_feat));

meg_feat0 = meg_feat(sig_feat_ind);
h_mat_meg0 = h_mat_meg(:,sig_feat_ind);

outpath_kmeans = fullfile(analysis_path,'kmeans');
if xlsout_save_flag
    if ~exist(outpath_kmeans,'dir'),mkdir(outpath_kmeans),end
end
outfile = fullfile(outpath_kmeans,'MEG_LGM_kmeans_sqeuclidean_k4.mat');
if ~exist(outfile,'file')
    nclust = 4;
    [clust_ind,clust_centroid,sumD,D] = clust_func(h_mat_meg0',nclust);
    
    if xlsout_save_flag
        save(outfile,'nclust','clust_ind','clust_centroid','sumD','D')
    end
else
    load(outfile)
    nclust = length(unique(clust_ind));
end

feat_group = {}; h_group = [];
c_group = {};
for iclust = 1:nclust
    feat_group0 = meg_feat0(clust_ind==iclust);
    feat_group{iclust,1} = feat_group0;
    
    h_group0 = h_mat_meg0(:,clust_ind==iclust);
    h_group(:,iclust) = mean(h_group0,2);
    
    mh_group0 = mh_mat_meg(clust_ind==iclust);
    std_h_group0 = sh_mat_meg(clust_ind==iclust);
    num_sigFDR_meg0 = num_sigFDR_meg(clust_ind==iclust);
    
    entry = compose('%0.2f (%0.2f)',mh_group0,std_h_group0);
    entry2 = compose('%0.2f',num_sigFDR_meg0);
    
    tok = cellfun(@(x) strsplit(x,'_'),feat_group0,'UniformOutput',0);
    c1 = vertcat(tok{:});
    c2 = [c1(:,[3,4,2]), entry, entry2];
    c_group{iclust,1} = c2;
end
R_fmri = corr(h_mat_fmri,clust_centroid');

% figure,imagesc(h_group,[0,1]),colorbar,colormap(jet)

mh_group = mean(h_group,1);
[~,hg_sort_ind] = sort(mh_group,'descend');
c_group_sort = c_group(hg_sort_ind);

c_tit = {'','Connectivity Metric','Frequency Band','Graph Measure','Mean (SD) Across Nodes','Prop. of Sig. Nodes'};
c_label = {};
for iclust = 1:nclust
    c1 = c_group_sort{iclust};
    c_label{iclust,1} = [{['Cluster ',num2str(iclust)]}; cell(size(c1,1)-1,1)];
end

cc_group = [c_tit; vertcat(c_label{:}),vertcat(c_group_sort{:})];

outfile_lgm = fullfile(outdir_fig,'LGM_tables.xlsx');
if xlsout_save_flag
    xlswrite(outfile_lgm,cc_group,'MEG_LGM_kmeans_sqeuclidean_k4')
end

%% make fMRI LGM h2 table

mh_mat_fmri = mean(h_mat_fmri,1)';
sh_mat_fmri = std(h_mat_fmri,[],1)';

entry_f = compose('%0.2f (%0.2f)',mh_mat_fmri,sh_mat_fmri);
entry2_f = compose('%0.2f',num_sigFDR_fmri);

tokf = cellfun(@(x) strsplit(x,'_'),fmri_feat,'UniformOutput',0);
cf1 = vertcat(tokf{:});
cf1_0 = replace(cf1,{'pos','neg','corr'},{'Positive','Negative','Cor'});
cf2 = [cf1_0(:,[3,4,2]), entry_f, entry2_f];

c_tit_f = {'Connectivity Metric','Positive/Negative Network','Graph Measure','Mean (SD) Across Nodes','Prop. of Sig. Nodes'};

c_group_f = [c_tit_f; cf2];

if xlsout_save_flag
    xlswrite(outfile_lgm,c_group_f,'fMRI_LGM')
end

%% make LGM h2 figures 
% 1 - MEG SOLAR, 2 - fMRI SOLAR
analysis_flags = 1;
for ia = 1:length(analysis_flags)
    analysis_flag = analysis_flags(ia);
    
    if analysis_flag == 1
        indir_root1 = fullfile(outpath_root,'SOLAR');
    elseif analysis_flag == 2
        indir_root1 = fullfile(outpath_root,'SOLAR_fMRI');
    end
    
    if ismember(analysis_flag,[1,3])
        % sel_conn = {'dwPLI','lcAEC'};
        sel_conn = {'dwPLI','AEC','lcAEC'};
        
        % sel_fnames = 'Alpha';
        sel_fnames = MEG_LGM_H2.freqBands;
        
        ftag = 'HCP_MEG_100%PropThres_';
        im_folder = 'bnet_png_FDR';
        
        tsp_opts = {length(sel_conn),4,0.01,0.02,0.02};
        
        outpath = fullfile(outpath_root,'SOLAR');
        outpath_LGM = fullfile(outpath,'_fig_a');
        if ~exist(outpath_LGM,'dir'),mkdir(outpath_LGM),end
    elseif ismember(analysis_flag,[2,4])
        varnames = replace(fmri_feat,'fMRI_','');
        ftag = 'HCP_fMRI_100%PropThres_';
        im_folder = 'bnet_png_FDR';
        
        tsp_opts = {2,4,0.01,0.02,0.02};
        
        outpath = fullfile(outpath_root,'SOLAR_fMRI');
        outpath_LGM = fullfile(outpath,'_fig_a');
        if ~exist(outpath_LGM,'dir'),mkdir(outpath_LGM),end
    end
    
    if analysis_flag == 1
        
        for ifr = 1:length(sel_fnames)
            sel_fname = sel_fnames{ifr};
            sel_var = contains(meg_feat,sel_conn) & contains(meg_feat,sel_fname);
            varnames = replace(meg_feat(sel_var),'MEG_','');
            
            fig = figure('Position',[50,1,1400,270*length(sel_conn)]);
            H = tight_subplot(tsp_opts{:});
            for ivar = 1:length(varnames)
                varname = varnames{ivar};
                tok = strsplit(varname,'_');
                graphMeasure = tok{1};
                tit0 = strjoin(strsplit(varname,'_'),' ');
                
                indir1 = fullfile(indir_root1,graphMeasure,im_folder);
                infile = fullfile(indir1,[ftag,varname,'_hvalGT_FDR.png']);
                im = imread(infile);
                im0 = padarray(im,[10,10],0);
                axes(H(ivar)),imshow(im0),axis off, axis image
                % text(750,725,tit0,'FontSize',8)
            end
            
            outfile = fullfile(outpath_LGM,[ftag,strjoin(sel_conn,'-'),'_',sel_fname,'_hvalGT_FDR']);
            dpi = '1200';
            img_outfile = outfile;
            img_outfile2 = [outfile,'_',dpi,'dpi'];
            save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
        end
        
    elseif analysis_flag == 2
        
        fig = figure('Position',[50,1,1400,540]);
        H = tight_subplot(tsp_opts{:});
        for ivar = 1:length(varnames)
            varname = varnames{ivar};
            tok = strsplit(varname,'_');
            graphMeasure = tok{1};
            tit0 = replace(varname,{'_','neg','pos','corr'},{' ','Negative','Positive','Cor'});
            
            indir1 = fullfile(indir_root1,graphMeasure,im_folder);
            infile = fullfile(indir1,[ftag,varname,'_hvalGT_FDR.png']);
            im = imread(infile);
            im0 = padarray(im,[10,10],0);
            axes(H(ivar)),imshow(im0),axis off, axis image
            % text(750,725,tit0,'FontSize',8)
        end
        
        outfile = fullfile(outpath_LGM,[ftag,'LGM','all','_hvalGT_FDR']);
        dpi = '1200';
        img_outfile = outfile;
        img_outfile2 = [outfile,'_',dpi,'dpi'];
        save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
    end
    
end

figure_name = ['colorbar3'];
fig = figure;
imagesc(checkerboard(10),[0,1])
colormap('jet')
colorbar
dpi = '1200';
img_outfile = [outdir_fig,'\',figure_name];
img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)

%% plot LGM h2 values across brain regions

[ix_list,ix_title] = atlas_Brainnetome_findLabelIndices_v2(nodelabel);
sel_ix_title = ix_list([4:9,3],1);
sel_ix_list = ix_list([4:9,3],[3,4]);
ix_limbic = {cell2mat(sel_ix_list([4,5],1)'),cell2mat(sel_ix_list([4,5],2)')};
sel_ix_title1 = {'Frontal'; 'Temporal'; 'Parietal'; 'Limbic'; 'Occipital'; 'Subcortical'};
sel_ix_list1 = [sel_ix_list(1:3,:); ix_limbic; sel_ix_list(6:7,:)];

sort_nodes0 = reshape(sel_ix_list1',[2*length(sel_ix_title1),1])';
sort_nodes = [sort_nodes0{:}]';

xtt = cellfun(@length,sort_nodes0,'UniformOutput',0)';
xTT = cumsum(cell2mat(xtt));
xTTT = [1 + [0; xTT(1:end-1)]; nnode];
xT4 = round((xTTT(1:end-1) + xTTT(2:end))/2);
xT5 = round(sum(reshape(xT4,[2,length(xT4)/2]),1)/2);

ix_title_lr0 = [strcat({'R.'},sel_ix_title1),strcat({'L.'},sel_ix_title1)]';
ix_title_lr = ix_title_lr0(:);

ix_tit_lr = repmat({'R','L'},[1,length(sel_ix_title1)]);

%%%%
outpath = fullfile(outpath_root,'SOLAR');
outpath1 = fullfile(outpath,'_fig3');
if ~exist(outpath1,'dir'),mkdir(outpath1),end

sel_conn = MEG_LGM_H2.conmethods([1,3]);
% sel_conn = MEG_LGM_H2.conmethods;

for ifr = 1:length(MEG_LGM_H2.freqBands)
    fname = MEG_LGM_H2.freqBands{ifr};
    
    h_mats = {}; tit_labs = {};
    for igraph = 1:length(MEG_LGM_H2.graphMeasures)
        gname = MEG_LGM_H2.graphMeasures{igraph};
        LGM_name = LGM_names{igraph};
        sel_ind = find(contains(meg_feat,gname) & contains(meg_feat,fname) & contains(meg_feat,sel_conn));
        h_mats{igraph,1} = h_mat_meg(sort_nodes,sel_ind);
        tit_labs{igraph} = [LGM_name,' (',gname,')'];
    end
    y_lab = 'Heritability (h^2)';
    
    fig = xx_lineplot(h_mats,length(MEG_LGM_H2.graphMeasures),xTTT,xT4,xT5,{},ix_tit_lr,sel_ix_title1,y_lab,tit_labs,sel_conn,0);
    
    if length(sel_conn) == 2
        outfile = fullfile(outpath1,['MEG_dwPLI-lcAEC_',fname,'_h2_across_ROI']);
    else
        outfile = fullfile(outpath1,['MEG_',fname,'_h2_across_ROI']);
    end
    dpi = '1200';
    img_outfile = outfile;
    img_outfile2 = [outfile,'_',dpi,'dpi'];
    save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
end

outpath = fullfile(outpath_root,'SOLAR_fMRI');
outpath1 = fullfile(outpath,'_fig3');
if ~exist(outpath1,'dir'),mkdir(outpath1),end

h_mats = {}; tit_labs = {};
for igraph = 1:length(fMRI_LGM_H2.graphMeasures)
    gname = fMRI_LGM_H2.graphMeasures{igraph};
    LGM_name = LGM_names{igraph};
    sel_ind = find(contains(fmri_feat,gname));
    h_mats{igraph,1} = h_mat_fmri(sort_nodes,sel_ind);
    tit_labs{igraph} = [LGM_name,' (',gname,')'];
end
y_lab = 'Heritability (h^2)';

fig = xx_lineplot(h_mats,length(fMRI_LGM_H2.graphMeasures),xTTT,xT4,xT5,{},ix_tit_lr,sel_ix_title1,y_lab,tit_labs,polarity_names,1);

outfile = fullfile(outpath1,['fMRI_','h2_across_ROI']);
dpi = '1200';
img_outfile = outfile;
img_outfile2 = [outfile,'_',dpi,'dpi'];
save_figure_hp(fig,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)

%% MEG LGM h2 (SOLAR vs APACE)

outpath_APACE = fullfile(outpath_root,'APACE');

pair_px_all = []; mdiff_all = []; rr_all = []; pair_ci_all = [];
tab_struc_meg = [];

mean_hstat_vals = []; % last row = [mean_val,low_ci,high_ci,pval]
vox_hvals = [];
vox_pvals_FWE = [];
vox_pvals_u = [];
vox_pvals_FDR = [];
nsig_all_FWE = [];
nsig_all_u = [];
nsig_all_FDR = [];
for igraph = 1:length(MEG_LGM_H2.graphMeasures)
    outpath0 = fullfile(outpath_APACE,MEG_LGM_H2.graphMeasures{igraph});
    for imethod = 1:length(MEG_LGM_H2.conmethods)
        for ifreq = 1:length(MEG_LGM_H2.freqBands)
            sel_measure = [MEG_LGM_H2.graphMeasures{igraph},'_',MEG_LGM_H2.conmethods{imethod},'_',MEG_LGM_H2.freqBands{ifreq}];
            outfolder = sel_measure;
            outpath_results = fullfile(outpath0,outfolder);
            
            p_file = fullfile(outpath_results,'AE_A_LRT_vox_FWEP.mat');
            pu_file = fullfile(outpath_results,'AE_A_LRT_vox_P.mat');
            pfdr_file = fullfile(outpath_results,'AE_A_LRT_vox_FDRP.mat');
            h_file = fullfile(outpath_results,'AE_A_h2.mat');
            pm_file = fullfile(outpath_results,'Pvals_h2.mat');
            ci_file = fullfile(outpath_results,'Boot_CIs.mat');
            
            % retrieve element-wise h2 and p values
            load(p_file)
            load(pu_file)
            load(pfdr_file)
            load(h_file)
            load(pm_file)
            load(ci_file)
            
            p_FWE = 10.^-AE_A_LRT_vox_FWEP;
            p_u = 10.^-AE_A_LRT_vox_P;
            p_FDR = 10.^-AE_A_LRT_vox_FDRP;
            
            nsig_FWE = length(find(p_FWE<0.05));
            nsig_u = length(find(p_u<0.05));
            nsig_FDR = length(find(p_FDR<0.05));
            
            mean_hstat_vals(igraph,imethod,ifreq,:) = [mean(AE_A_h2),CIs_h2(1,1),CIs_h2(1,2),Pvals_h2(1)];
            vox_hvals(igraph,imethod,ifreq,:) = AE_A_h2;
            vox_pvals_FWE(igraph,imethod,ifreq,:) = p_FWE;
            vox_pvals_u(igraph,imethod,ifreq,:) = p_u;
            vox_pvals_FDR(igraph,imethod,ifreq,:) = p_FDR;
            nsig_all_FWE(igraph,imethod,ifreq) = nsig_FWE;
            nsig_all_u(igraph,imethod,ifreq) = nsig_u;
            nsig_all_FDR(igraph,imethod,ifreq) = nsig_FDR;
            
            ind = contains(MEG_LGM_H2.tit_lgm,sel_measure);
            sel_tit = MEG_LGM_H2.tit_lgm(ind)';
            h_solar = MEG_LGM_H2.h_all(ind);
            [hx,px,ci_x] = ttest(h_solar,AE_A_h2);
            pair_px_all(igraph,imethod,ifreq) = px;
            pair_ci_all(igraph,imethod,ifreq,:) = ci_x;
            
            mdiff_all(igraph,imethod,ifreq,:) = AE_A_h2-h_solar;
            [rx,prx] = corr(h_solar,AE_A_h2);
            rr_all(igraph,imethod,ifreq) = rx;
            
            entry_d = compose('r = %.2f, d = [%.2f,%.2f]',rx,ci_x(1),ci_x(2));
            tab_struc_meg{imethod,1}{igraph,ifreq} = entry_d{1};
        end
    end
end

tab_struc_meg0 = {};
for imethod = 1:length(MEG_LGM_H2.conmethods)
    blank_tit = cell(1,nfreq+1);
    blank_tit{1} = MEG_LGM_H2.conmethods{imethod};
    tab_struc_meg0{imethod,1} = [blank_tit; [{[]},MEG_LGM_H2.freqBands]; [MEG_LGM_H2.graphMeasures',tab_struc_meg{imethod}]];
end
tab_struc_meg0 = vertcat(tab_struc_meg0{:});

outfile_s1 = fullfile(outdir_fig,'LGM_SOLAR-vs-APACE.xlsx');
if xlsout_save_flag
    xlswrite(outfile_s1,tab_struc_meg0,'MEG')
end

%% fMRI LGM h2 (SOLAR vs APACE)

outpath_APACE = fullfile(outpath_root,'APACE_fMRI');

pair_px_all = []; mdiff_all = []; rr_all = [];
tab_struc_fmri = [];

mean_hstat_vals = []; % last row = [mean_val,low_ci,high_ci,pval]
vox_hvals = [];
vox_pvals_FWE = [];
vox_pvals_u = [];
vox_pvals_FDR = [];
nsig_all_FWE = [];
nsig_all_u = [];
nsig_all_FDR = [];
for igraph = 1:length(fMRI_LGM_H2.graphMeasures)
    outpath0 = fullfile(outpath_APACE,fMRI_LGM_H2.graphMeasures{igraph});
    for ipol = 1:length(fMRI_LGM_H2.polarity_all)
        sel_measure = [fMRI_LGM_H2.graphMeasures{igraph},'_',fMRI_LGM_H2.conmethods{1},'_',fMRI_LGM_H2.polarity_all{ipol}];
        
        ind = contains(fMRI_LGM_H2.tit_lgm,sel_measure);
        sel_tit = fMRI_LGM_H2.tit_lgm(ind)';
        
        outfolder = sel_measure;
        outpath_results = fullfile(outpath0,outfolder);
        
        p_file = fullfile(outpath_results,'AE_A_LRT_vox_FWEP.mat');
        pu_file = fullfile(outpath_results,'AE_A_LRT_vox_P.mat');
        pfdr_file = fullfile(outpath_results,'AE_A_LRT_vox_FDRP.mat');
        h_file = fullfile(outpath_results,'AE_A_h2.mat');
        pm_file = fullfile(outpath_results,'Pvals_h2.mat');
        ci_file = fullfile(outpath_results,'Boot_CIs.mat');
        
        % retrieve element-wise h2 and p values
        load(p_file)
        load(pu_file)
        load(pfdr_file)
        load(h_file)
        load(pm_file)
        load(ci_file)
        
        p_FWE = 10.^-AE_A_LRT_vox_FWEP;
        p_u = 10.^-AE_A_LRT_vox_P;
        p_FDR = 10.^-AE_A_LRT_vox_FDRP;
        
        nsig_FWE = length(find(p_FWE<0.05));
        nsig_u = length(find(p_u<0.05));
        nsig_FDR = length(find(p_FDR<0.05));
        
        mean_hstat_vals(igraph,ipol,:) = [mean(AE_A_h2),CIs_h2(1,1),CIs_h2(1,2),Pvals_h2(1)];
        vox_hvals(igraph,ipol,:) = AE_A_h2;
        
        vox_pvals_FWE(igraph,ipol,:) = p_FWE;
        vox_pvals_u(igraph,ipol,:) = p_u;
        vox_pvals_FDR(igraph,ipol,:) = p_FDR;
        
        nsig_all_FWE(igraph,ipol) = nsig_FWE;
        nsig_all_u(igraph,ipol) = nsig_u;
        nsig_all_FDR(igraph,ipol) = nsig_FDR;
        
        h_solar = fMRI_LGM_H2.h_all(ind);
        [hx,px,ci_x] = ttest(h_solar,AE_A_h2);
        pair_px_all(igraph,ipol) = px;
        
        mdiff_all(igraph,ipol,:) = AE_A_h2-h_solar;
        [rx,prx] = corr(h_solar,AE_A_h2);
        rr_all(igraph,ipol) = rx;
        
        entry_d = compose('r = %.2f, d = [%.2f,%.2f]',rx,ci_x(1),ci_x(2));
        tab_struc_fmri{igraph,ipol} = entry_d{1};
    end
    
end

tab_struc_fmri0 = [[{[]},polarity_names]; [fMRI_LGM_H2.graphMeasures',tab_struc_fmri]];

if xlsout_save_flag
    xlswrite(outfile_s1,tab_struc_fmri0,'fMRI')
end

%% load meg GGM h2 values (sensor space)
MEG_GGM_H2_Source = MEG_GGM_H2;

indir_meg_ggm_sensor = fullfile(indir,'phaatef_MEG_GGM_Sensor_genetictest');
outfile = fullfile(outpath_root,'SOLAR','MEG_GGM_Sensor_100%PropThres_hval.mat');

if ~exist(outfile,'file')
    infile_GGM = fullfile(indir_meg_ggm_sensor,'HCP_MEG_GGM_Sensor_100%PropThres.csv');
    tab = readtable(infile_GGM);
    
    tit = tab.Properties.VariableNames;
    tit_ggm = tit(7:end);
    
    conmethods = {'dwPLI','AEC','lcAEC'};
    freqBands = {'Delta','Theta','Alpha','lBeta','hBeta','lGamma'};
    graphMeasures = {'GE','CPL','T','S'};
    
    h_all = []; p_all = []; he_all = []; cov_all = []; kurt_all = []; kurt_score_all = {};
    for iggm = 1:length(tit_ggm)
        varname = tit_ggm{iggm};
        indir0 = fullfile(indir_meg_ggm_sensor,varname);
        
        [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0);
        
        h_all(iggm,1) = hval;
        p_all(iggm,1) = pval;
        he_all(iggm,1) = herror;
        cov_all(iggm,1) = covar_var;
        kurt_all(iggm,1) = kurt;
        kurt_score_all{iggm,1} = kurt_score;
    end
    
    pFDR_all = mafdr(p_all,'BHFDR',1);
    
    MEG_GGM_H2 = [];
    MEG_GGM_H2.h_all = h_all;
    MEG_GGM_H2.p_all = p_all;
    MEG_GGM_H2.pFDR_all = pFDR_all;
    MEG_GGM_H2.he_all = he_all;
    MEG_GGM_H2.cov_all = cov_all;
    MEG_GGM_H2.kurt_all = kurt_all;
    MEG_GGM_H2.kurt_score_all = kurt_score_all;
    MEG_GGM_H2.tit_ggm = tit_ggm;
    MEG_GGM_H2.conmethods = conmethods;
    MEG_GGM_H2.freqBands = freqBands;
    MEG_GGM_H2.graphMeasures = graphMeasures;
    
    save(outfile,'MEG_GGM_H2')
else
    load(outfile)
end

MEG_GGM_H2_Sensor = MEG_GGM_H2;
clear MEG_GGM_H2

nmethod = length(MEG_GGM_H2_Source.conmethods);
nfreq = length(MEG_GGM_H2_Source.freqBands);
ngraph = length(MEG_GGM_H2_Source.graphMeasures);
GGM_names = {'Global Efficiency','Characteristic Path Length','Transitivity','Synchronizability'};

% plot supplementary figure (h2 sensor vs source space)
figPos_b = [10,100,screenPos(3)*0.975,screenPos(4)*0.6];

space_tit = {'Source','Sensor'};
H2 = []; 
H2(1,:,:,:) = reshape(MEG_GGM_H2_Source.h_all,[ngraph,nmethod,nfreq]);
H2(2,:,:,:) = reshape(MEG_GGM_H2_Sensor.h_all,[ngraph,nmethod,nfreq]);
theta = deg2rad(0:360/6:360);

for imethod = 1:size(H2,3)
    conmethod_i = MEG_GGM_H2_Source.conmethods{imethod};
    file_name = ['H2_SensorVsSource_',conmethod_i];
    figure_name = ['H2','_',conmethod_i];
    h = figure('Position',figPos_b);
    set(h,'name',figure_name)
    st = suptitle(['Heritability of Global Graph Measures for MEG ',conmethod_i]);
    set(st,'FontSize',16)
    
    for igraph = 1:size(H2,2)
        H2i = squeeze(H2(:,igraph,imethod,:));
        subplot(1,4,igraph);
        polarplot(theta,[H2i,H2i(:,1)],'.-','LineWidth',1.5,'MarkerSize',14);
        rlim([0,1])
        thetaticks(0:360/6:360-360/6)
        thetaticklabels(MEG_GGM_H2_Source.freqBands)
        title(GGM_names{igraph})
        
        pos_s = get(gca,'Position');
        if igraph == 1, newpos_sx = pos_s(1)-0.09;
        elseif igraph == 2, newpos_sx = pos_s(1)-0.05; middle_pos = newpos_sx;
        elseif igraph == 3, newpos_sx = pos_s(1)-0.01; middle_pos = (middle_pos + newpos_sx)/2;
        elseif igraph == 4, newpos_sx = pos_s(1)+0.03; end
        
        newpos_sy = pos_s(2)-0.04;
        newpos_s = [newpos_sx newpos_sy pos_s(3)*1.1 pos_s(4)*1.1];
        set(gca,'Position',newpos_s,'FontSize',13,'GridAlpha',0.2,'GridColor','k')
    end
    
    lgd = legend(strcat(space_tit,' Space'),'Orientation','Horizontal','Location','southoutside');
    pos_l = get(lgd,'Position');
    set(lgd,'Position',[middle_pos pos_l(2)*0.5 pos_l(3) pos_l(4)],'FontSize',14)
    
    img_outfile = [outdir_fig,'\',file_name];
    dpi = '1200';
    img_outfile2 = [outdir_fig,'\',file_name,'_',dpi,'dpi'];
    save_figure_hp(h,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
end

