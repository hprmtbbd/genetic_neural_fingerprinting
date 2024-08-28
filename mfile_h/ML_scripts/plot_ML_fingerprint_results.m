
clear,clc

% written for MATLAB R2023b (needed newer version for swarmchat and boxchart)

workspace_path = 'D:\workspace';
proj_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript');
analysis_path = fullfile(proj_path,'h_analysis');

func_path = fullfile(workspace_path,'MATLAB','functions');
addpath(genpath(func_path))

outdir_root = 'D:\Data\h_analysis_out';
outdir_fig = fullfile(proj_path,'Figures','version3_revision1');
if ~exist(outdir_fig,'dir'),mkdir(outdir_fig),end

outfolder = 'ML_fingerprint_v2';
outfolder2 = 'featsel_0';

%%
col_tit = {'TPR','TNR','PPV','AUC'};
tit_a = {'MZ','SIB','UR'};
tits = {'MZ twins','Non-MZ siblings','Unrelated subjects','Weighted average'}';

screenPos = get(0,'Screensize');
save_flag_png = 0; save_flag_fig = 0;

xlsout_flag = 0;

data_tags = {'meg','fmri'};
measure_tags = {'GGM','LGM'};
measure_tits = {'Global Graph Measures','Local Graph Measures'};

%%
analysis_flags = 1:4; % 1:4 [1,2]

acc_array = []; auc_array = [];
acc_mat = []; auc_mat = [];
acc_s_mat = []; auc_s_mat = [];
acc_p_mat = []; auc_p_mat = [];
class_tits = {}; class_tags = {};
for ia = 1:length(analysis_flags)
    analysis_flag = analysis_flags(ia);
    
    if analysis_flag == 1
        class_tag = 'threeClass';
        class_tit = 'Three group classification';
    elseif analysis_flag == 2
        class_tag = [tit_a{1},'-',tit_a{3}];
        class_tit = [tits{1},' vs. ',tits{3}];
    elseif analysis_flag == 3
        class_tag = [tit_a{1},'-',tit_a{2}];
        class_tit = [tits{1},' vs. ',tits{2}];
    elseif analysis_flag == 4
        class_tag = [tit_a{2},'-',tit_a{3}];
        class_tit = [tits{2},' vs. ',tits{3}];
    end
    
    for imeasure = 1:length(measure_tags)
        measure_tag0 = measure_tags{imeasure};
        
        if strcmp(measure_tag0,'GGM')
            model_tag = 'SVM-linear-1-abs';
        elseif strcmp(measure_tag0,'LGM')
            model_tag = 'SVM-linear-1-corr';
        end
        
        acc_array0 = []; auc_array0 = [];
        acc_vec = []; auc_vec = [];
        acc_s_vec = []; auc_s_vec = [];
        acc_p_vec = []; auc_p_vec = [];
        
        for idata = 1:length(data_tags)
            data_tag0 = data_tags{idata};
            
            outdir = fullfile(outdir_root,outfolder,outfolder2,[measure_tag0,'_',data_tag0]);
            outfile = fullfile(outdir,[model_tag,'_',class_tag,'.mat']);
            outfile_perm = fullfile(outdir,[model_tag,'_',class_tag,'_perm1000.mat']);
            load(outfile)
            load(outfile_perm)
            
            nmethod = length(ML_STATS_all);
            for imethod = 1:nmethod
                ML_STATS = ML_STATS_all{imethod};
                
                acc_m = squeeze(nanmean(ML_STATS.acc_all,[1,2]));
                auc_m = squeeze(nanmean(ML_STATS.auc_all(:,:,end),[1,2]));
                acc_s = squeeze(nanstd(ML_STATS.acc_all,[],[1,2]));
                auc_s = squeeze(nanstd(ML_STATS.auc_all(:,:,end),[],[1,2]));
                
                m_acc_perm = squeeze(mean(ML_STATS_perm{imethod}.acc_perm,[1,2]));
                m_auc_perm = squeeze(mean(ML_STATS_perm{imethod}.auc_perm(:,:,end,:),[1,2]));
                
                [sig_acc,px_acc,pj_acc,xj_acc] = perm_test(acc_m,m_acc_perm,0.05,1);
                [sig_auc,px_auc,pj_auc,xj_auc] = perm_test(auc_m,m_auc_perm,0.05,1);
                
                acc_vec = [acc_vec; acc_m];
                auc_vec = [auc_vec; auc_m];
                
                acc_s_vec = [acc_s_vec; acc_s];
                auc_s_vec = [auc_s_vec; auc_s];
                
                acc_p_vec = [acc_p_vec; px_acc];
                auc_p_vec = [auc_p_vec; px_auc];
                
                acc_array0 = cat(3,acc_array0,ML_STATS.acc_all);
                auc_array0 = cat(3,auc_array0,ML_STATS.auc_all(:,:,end));
            end
        end
        
        acc_array(ia,imeasure,:,:,:) = permute(acc_array0,[3,1,2]); %(iclass,imeasure,imethod,iperm,ifold)
        auc_array(ia,imeasure,:,:,:) = permute(auc_array0,[3,1,2]);
        
        acc_mat(ia,imeasure,:) = acc_vec; auc_mat(ia,imeasure,:) = auc_vec;
        acc_s_mat(ia,imeasure,:) = acc_s_vec; auc_s_mat(ia,imeasure,:) = auc_s_vec;
        acc_p_mat(ia,imeasure,:) = acc_p_vec; auc_p_mat(ia,imeasure,:) = auc_p_vec;
    end
    
    class_tits{ia,1} = class_tit; class_tags{ia,1} = class_tag;
end

method_tits = {'MEG dwPLI'; 'MEG AEC'; 'MEG lcAEC'; 'fMRI Cor'};
metric_tits = {'Accuracy','AUC'};

%%
% figPos_a = [50,20,screenPos(3)*0.2*length(measure_tits),screenPos(4)];
figPos_a = [50,50,screenPos(3)*0.8,screenPos(4)*0.5];

alpha = [0.05,0.005];

for iclass = 1:length(class_tits)
    class_tag0 = class_tags{iclass};
    class_tit0 = class_tits{iclass};
    
    acc_mat0 = squeeze(acc_mat(iclass,:,:));
    auc_mat0 = squeeze(auc_mat(iclass,:,:));
    acc_s_mat0 = squeeze(acc_s_mat(iclass,:,:));
    auc_s_mat0 = squeeze(auc_s_mat(iclass,:,:));
    acc_p_mat0 = squeeze(acc_p_mat(iclass,:,:));
    auc_p_mat0 = squeeze(auc_p_mat(iclass,:,:));
    
    m_mat = cat(3,acc_mat0,auc_mat0);
    se_mat = cat(3,acc_s_mat0,auc_s_mat0);
    p_mat = cat(3,acc_p_mat0,auc_p_mat0);
    
    nmeasure = size(m_mat,1);
    nmethod = size(m_mat,2);
    nmetric = size(m_mat,3);
    
    % (imeasure,imetric,imethod)
    m = permute(m_mat,[1,3,2]);
    se = permute(se_mat,[1,3,2]);
    p = permute(p_mat,[1,3,2]);
    
    %-------------------------------------------------------
    figure_name = [class_tag0,'_ML'];
    
    h = figure;
    set(h,'name',figure_name,'position',figPos_a)
    for imetric = 1:nmetric
        fs = subplot(1,2,imetric);
        
        mu = squeeze(m(:,imetric,:));
        sigma = squeeze(se(:,imetric,:));
        p_values = squeeze(p(:,imetric,:));
        
        p_values_adj = nan([size(p_values),nmethod]);
        for imethod = 1:nmethod
            p_values_adj(:,imethod,imethod) = p_values(:,imethod);
        end
        
        p_temp = p_values_adj;
        N1 = size(p_temp,1);
        N2 = size(p_temp,2);
        
        p_val = [];
        for i = 1:N2
            p_val(:,i) = p_temp(:,i,i);
        end
        
        C = [0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250;  0.4660 0.6740 0.1880; 0 0.4470 0.7410];
        
        ax = superbar(mu,'BarFaceColor',permute(C,[3,1,2]),'E',sigma,'ErrorbarColor','k','P',p_val,...
        'PStarShowGT',false,'PStarThreshold',alpha,'PStarIcon','*',...
        'PStarShowNS',false,'PStarFontSize',14,'PLineColor','k','PLineWidth',2,'PStarColor','k',...
        'PStarOffset',0.02);
        
        [m1,n1] = size(mu');
        loc1 = 0.15* linspace(-ceil(m1/2),ceil(m1/2),m1)'*ones(1,n1) + ones(m1,1)*(1:n1);
        xlim( [ 0.5 ( max(loc1(:))+min(loc1(:))/2 ) ] )
        xticks(1:nmetric)
        set(gca,'XTickLabel',measure_tits,'FontName','Arial')
        set(gca,'Fontsize',12)
        ylim([0,1.09]),yticks(0:0.1:1)
        % ylabel(metric_tits{imetric})
        set(gca,'YGrid','on')
        tsub = title(metric_tits{imetric});
        set(tsub,'FontSize',14)
        
        pos_s = get(fs,'Position');
        
        newpos_s = [pos_s(1) pos_s(2)+0.03 pos_s(3) pos_s(4)*0.9];
        set(fs,'Position',newpos_s)
    end
    
    lgd = legend(ax(1,:),method_tits,'Orientation','Horizontal');
    pos_l = get(lgd,'Position');
    set(lgd,'Position',[pos_l(1)-0.18 pos_l(2)*0.2-0.14 pos_l(3) pos_l(4)])
    
    sgt = annotation('textbox','String',class_tit0,'EdgeColor','none','FontSize',18);
    set(sgt,'FitBoxToText','on')
    drawnow
    tl_pos = get(sgt,'Position'); xtl = 0.5-tl_pos(3)/2+0.035;
    set(sgt,'Position',[xtl 1 1 0])
    
    dpi = '1200';
    img_outfile = [outdir_fig,'\',figure_name];
    img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(h,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)
end

%% figures (second version)

figPos_b = [50,50,screenPos(3)*0.8,screenPos(4)*0.5];
C = [0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250;  0.4660 0.6740 0.1880; 0 0.4470 0.7410];

for iclass = 2 % 1:length(class_tits)
    class_tag0 = class_tags{iclass};
    figure_name = [class_tag0,'_ML_violinplot'];

    h = figure('Position',figPos_b);
    T = tiledlayout(1,2);
    for imetric = 1:2
        t = tiledlayout(T,1,length(measure_tits));
        t.Layout.Tile = imetric;

        if imetric == 1
            metric_array = acc_array;
            metric_p_mat = acc_p_mat;
            metric_tit = 'Accuracy';
        elseif imetric == 2
            metric_array = auc_array;
            metric_p_mat = auc_p_mat;
            metric_tit = 'AUC';
        end

        for imeasure = 1:length(measure_tits)
            % subplot(1,length(measure_tits),imeasure)

            nexttile(t,imeasure)
            metric_array1 = squeeze(mean(metric_array(iclass,imeasure,:,:,:),5));
            x_vec = 1.25*(1:4)';
            xx_vec = mat2cell([x_vec-0.3,x_vec+0.3],ones(size(x_vec)),2);
            xx = x_vec.*ones(size(metric_array1));
            metric_p1 = squeeze(metric_p_mat(iclass,imeasure,:));

            for imethod = 1:length(method_tits)
                swarmchart(xx(imethod,:),metric_array1(imethod,:),10,C(imethod,:),'filled',...
                    'XJitterWidth',0.8*(min(diff(x_vec))))
                hold on
            end
            xlim([0.5,5.75]),ylim([0.4,1.06])
            xticks([])
            xlabel(measure_tits{imeasure})
            ps = sigstar(xx_vec,metric_p1);

            set(gca,'FontSize',12)
            set(ps(:,2),'FontSize',14)

            % legend(method_tits)
        end
        t.TileSpacing = 'compact';
        t.Title.String = [metric_tit, ' (5-fold Cross-Validation)'];
        t.Title.FontWeight = 'bold';
        t.Title.FontSize = 14;
    end
    T.TileSpacing = 'compact';
    T.Title.String = class_tits{iclass};
    T.Title.FontSize = 18;

    ax = axes(T,'Visible','off');
    for imethod = 1:length(method_tits)
        swarmchart(ax,NaN,NaN,10,C(imethod,:),'filled')
        hold on
    end
    hold off
    set(ax,'Visible','off')

    lgd = legend(method_tits,'Orientation','horizontal','FontSize',12);
    lgd.Layout.Tile = 'south';

    dpi = '1200';
    img_outfile = [outdir_fig,'\',figure_name];
    img_outfile2 = [outdir_fig,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(h,save_flag_png,img_outfile2,dpi,save_flag_fig,img_outfile)

end

%% supp info

supp_dir = 'D:\Data\h_analysis_out\Supp_Data';
supp_dir1 = fullfile(supp_dir,'data4');
if ~exist(supp_dir1,'dir'),mkdir(supp_dir1),end

method_tits0 = replace(method_tits,' ','_');

for iclass = 1:length(class_tags)
    for imeasure = 1:length(measure_tits)
        ACC = []; AUC = [];
        for imethod = 1:length(method_tits)
            method_tit0 = method_tits0{imethod};
            ACC.(method_tit0) = squeeze(acc_array(iclass,imeasure,imethod,:,:));
            AUC.(method_tit0) = squeeze(auc_array(iclass,imeasure,imethod,:,:));
        end
        ACC.p_value = array2table(squeeze(acc_p_mat(iclass,imeasure,:))','VariableNames',method_tits0');
        AUC.p_value = array2table(squeeze(auc_p_mat(iclass,imeasure,:))','VariableNames',method_tits0');

        supp_file = fullfile(supp_dir1,[measure_tags{imeasure},'_',class_tags{iclass},'.mat']);
        if ~exist(supp_file,'file')
            save(supp_file,'ACC','AUC')
        end
    end
end

%% tables
outfile_xls = fullfile(outdir_fig,['ML_tables','.xlsx']);

for ia = 1:length(analysis_flags)
    analysis_flag = analysis_flags(ia);
    
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
    
    cc_all = {};
    for imeasure = 1:length(measure_tags)
        measure_tag0 = measure_tags{imeasure};
        measure_tit0 = measure_tits{imeasure};
        
        if strcmp(measure_tag0,'GGM')
            model_tag = 'SVM-linear-1-abs';
        elseif strcmp(measure_tag0,'LGM')
            model_tag = 'SVM-linear-1-corr';
        end
        
        c_all = {}; row_titc = {};
        for idata = 1:length(data_tags)
            data_tag0 = data_tags{idata};
            
            outdir = fullfile(outdir_root,outfolder,outfolder2,[measure_tag0,'_',data_tag0]);
            outfile = fullfile(outdir,[model_tag,'_',class_tag,'.mat']);
            load(outfile)
            nmethod = length(ML_STATS_all);
            
            for imethod = 1:nmethod
                ML_STATS = ML_STATS_all{imethod};
                
                conn_tit = method_tits{3*(idata-1)+imethod};
                con_tits = cell(1,4);
                con_tits(1) = compose('%s (ACC = %0.2f +/ %0.2f)',conn_tit,squeeze(nanmean(ML_STATS.acc_all,[1,2])),squeeze(nanstd(ML_STATS.acc_all,[],[1,2])));
                
                c1 = compose('%0.2f +/- %0.2f',squeeze(nanmean(ML_STATS.tpr_all,[1,2])),squeeze(nanstd(ML_STATS.tpr_all,[],[1,2])));
                c2 = compose('%0.2f +/- %0.2f',squeeze(nanmean(ML_STATS.tnr_all,[1,2])),squeeze(nanstd(ML_STATS.tnr_all,[],[1,2])));
                c3 = compose('%0.2f +/- %0.2f',squeeze(nanmean(ML_STATS.ppv_all,[1,2])),squeeze(nanstd(ML_STATS.ppv_all,[],[1,2])));
                c4 = compose('%0.2f +/- %0.2f',squeeze(nanmean(ML_STATS.auc_all,[1,2])),squeeze(nanstd(ML_STATS.auc_all,[],[1,2])));
                cc = [con_tits; col_tit; c1,c2,c3,c4];
                
                c_all = [c_all; cc];
                row_titc = [row_titc; [{[]; []}; row_tit]];
            end
        end
        
        measure_titc = cell(1,4); measure_titc(1) = {measure_tit0};
        cc_all = [cc_all,[measure_titc; c_all]];
    end
    
    ccc_all = [[{[]}; row_titc],cc_all];
    
    if xlsout_flag
        writecell(ccc_all,outfile_xls,'Sheet',class_tag)
    end
end
