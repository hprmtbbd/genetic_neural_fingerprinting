
clear,clc

% written for MATLAB R2023b (needed newer version for swarmchat and boxchart)

workspace_path = 'D:\workspace';
proj_path = fullfile(workspace_path,'MEGConnHeritability_Manuscript');
analysis_path = fullfile(proj_path,'h_analysis');

toolbox_path = fullfile(workspace_path,'MATLAB','toolbox');
func_path = fullfile(workspace_path,'MATLAB','functions');

excel_path = analysis_path;

addpath('D:\workspace\MEGConnHeritability_Manuscript\mfile_h\1',...
    genpath(func_path))

outdir = 'D:\workspace\MEGConnHeritability_Manuscript\Figures\version3\Analysis Pipeline';
if ~exist(outdir,'dir'),mkdir(outdir),end

%%
save_flag = 1;
data_flags = [1,2]; % 1 - meg, 2 - fmri
measure_flags = [1,2]; % 1 - global, 2 - local

%%
VV_mat = {};
pp_mat = [];

for imeasure = 1:length(measure_flags)
    measure_flag = measure_flags(imeasure);

    % select MEG or fMRI
    if measure_flag == 1
        measure_tag = 'GGM';
        measure_tit = 'Global Graph Measure';
        Gnames = {'GE','CPL','T','S'};
    elseif measure_flag == 2
        measure_tag = 'LGM';
        measure_tit = 'Local Graph Measure';
        Gnames = {'STR','EVC','CC','NE'};
        nnode = 246;
    end

    V_mat = {}; fig_tits = {};
    for idata = 1:length(data_flags)
        data_flag = data_flags(idata);

        if data_flag == 1
            infile = fullfile(analysis_path,['HCP_MEG_',measure_tag,'_Source_100%PropThres.csv']);
            data_tag = 'MEG';

            connames = {'dwPLI','AEC','lcAEC'};
            con_tits = connames;
            freqBands = {'Delta','Theta','Alpha','lBeta','hBeta','lGamma'};
            nfreq = length(freqBands);
        elseif data_flag == 2
            infile = fullfile(analysis_path,['HCP_fMRI_',measure_tag,'_Source_100%PropThres.csv']);
            data_tag = 'fMRI';

            polarity_all = {'pos','neg'};
            connames = {'corr'};
            con_tits = {'Cor'};
            npol = length(polarity_all);
        end

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
        ii_mz = [jmz, imz]';
        jj_mz = ii_mz(:);

        jj_single = find(all(phi0==0,2));

        if data_flag == 1
            [isib,jsib] = find(phi0_liu==1/2);

            ii_sib = [jsib, isib]';
            jj_sib = ii_sib(:);

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

            k_sib = cellfun(@(x) nchoosek(x,2),fam_locs00,'UniformOutput',0);
            jj_sib = cell2mat(k_sib)';
            isib = jj_sib(1,:)';
            jsib = jj_sib(2,:)';

            % isib = cell2mat(cellfun(@(x) x(1),fam_locs00,'UniformOutput',0));
            % jsib = cell2mat(cellfun(@(x) x(2),fam_locs00,'UniformOutput',0));

            nmz = length(jmz); nsib = length(jsib);
        end

        %%
        if 0 % data_flag == 1
            if measure_flag == 1
                sel_feat = [3,1,1];
                dat = m_ggm(:,sel_feat(1),sel_feat(2),sel_feat(3));
            elseif measure_flag == 2
                sel_feat = [3,4,1];
                dat = m_ggm(:,:,sel_feat(1),sel_feat(2),sel_feat(3));
            end
            tit = m_tit(sel_feat(1),sel_feat(2),sel_feat(3));
            pdat = (dat - mean(dat))./std(dat);

            jj = [jj_mz; jj_sib; jj_single];

            if measure_flag == 1
                pdat_sort = pdat(jj);
                D_sort = abs(pdat_sort-pdat_sort');
            elseif measure_flag == 2
                pdat_sort = pdat(jj,:);
                D_sort = pdist2(pdat_sort,pdat_sort,'correlation');
                R_sort = corr(pdat_sort,pdat_sort);
            end

            sel_inds = [21:22,7:8,13:14];
            R_sel = D_sort(sel_inds,sel_inds);
            R_sel0 = (R_sel-min(R_sel(:)))./(max(R_sel(:))-min(R_sel(:)));
            Itriu = ones(length(R_sel)); % triu(ones(length(R_sel)));
            R_sel0(~Itriu) = NaN;

            % labs = strcat({'MZ Twin Pair '},{'1a','1b','2a','2b','3a','3b'});
            labs = strcat({'MZ '},{'1a','1b','2a','2b','3a','3b'});

            figure_name = 'example_Graph_measure_distance_matrix';

            screenPos = get(0,'ScreenSize');
            fig = figure('Position',[20,20,screenPos(3)*0.5,screenPos(4)*0.8]);
            h = heatmap(labs,labs,R_sel0,'CellLabelColor','none');
            h.ColorLimits = [min(R_sel0(:)),max(R_sel0(:))];
            h.FontSize = 20;
            h.ColorbarVisible = 'off';
            s = struct(h);
            s.XAxis.TickLabelRotation = 0;
            % title('\fontsize{20} \bf{Input features: Graph measure distance scores}');
            colormap(brewermap([],'Reds'))

            img_outfile = fullfile(outdir,figure_name);
            if save_flag
                print(fig,img_outfile,'-dpng','-r600')
                saveas(fig,img_outfile)
            end
        end

        %%
        pfeats = (m_ggm - mean(m_ggm,1))./std(m_ggm,[],1);

        if measure_flag == 1
            pfeats_sort = pfeats(:,:);
            Dmat = [];
            for ifeat = 1:size(pfeats_sort,2)
                D_sort = abs(pfeats_sort(:,ifeat)-pfeats_sort(:,ifeat)');
                Dmat(:,:,ifeat) = D_sort;
            end
        elseif measure_flag == 2
            pfeats_sort = pfeats(:,:,:);
            Dmat = []; Rmat = [];
            for ifeat = 1:size(pfeats_sort,3)
                pfeat = pfeats_sort(:,:,ifeat);
                D_sort = pdist2(pfeat,pfeat,'correlation');
                R_sort = corr(pfeat',pfeat');
                Dmat(:,:,ifeat) = D_sort;
                Rmat(:,:,ifeat) = R_sort;
            end
        end

        mtit_vec = m_tit(:);

        vmat = [];
        for imethod = 1:length(connames)
            conmethod = connames{imethod};
            sel_ind = contains(mtit_vec,['_',conmethod,'_']);
            sel_Dmat = Dmat(:,:,sel_ind);
            vmat0 = sqrt(sum(sel_Dmat.^2,3));
            vmat(:,:,imethod) = vmat0;
        end

        sel_vmat = permute(vmat,[3,1,2]);

        all_lind = (1:nsubj*nsubj)';
        self_lind = (1:(nsubj+1):nsubj*nsubj)';
        mz_lind = sub2ind([nsubj,nsubj],imz,jmz);
        sib_lind = sub2ind([nsubj,nsubj],isib,jsib);
        ur_lind = setdiff(all_lind,[self_lind; mz_lind; sib_lind]);

        nur = length(ur_lind);

        mz_Dmat = sel_vmat(:,mz_lind);
        sib_Dmat = sel_vmat(:,sib_lind);
        ur_Dmat = sel_vmat(:,ur_lind);

        for imethod = 1:length(connames)
            V_mat = [V_mat; [{mz_Dmat(imethod,:)},{sib_Dmat(imethod,:)},{ur_Dmat(imethod,:)}]];
            fig_tits = [fig_tits; [data_tag,' ',con_tits{imethod}]];
        end
    end

    %%
    
    p_mat = [];
    for isel_feat = 1:size(V_mat,1)
        mz_Dmat0 = V_mat{isel_feat,1};
        sib_Dmat0 = V_mat{isel_feat,2};
        ur_Dmat0 = V_mat{isel_feat,3};
        p12 = ranksum(mz_Dmat0,sib_Dmat0);
        p13 = ranksum(mz_Dmat0,ur_Dmat0);
        p23 = ranksum(sib_Dmat0,ur_Dmat0);
        pvec = [p12,p13,p23];
        p_mat(isel_feat,:) = pvec;
    end

    VV_mat(:,:,imeasure) = V_mat;
    pp_mat(:,:,imeasure) = p_mat;
end

%%
group_labels = {'MZ Twins',['Non-MZ\newline',char(8201),'Siblings'],'Unrelated\newline Subjects'};
groups = {[1,2],[1,3],[2,3]};

pm_vec = reshape(pp_mat,[size(VV_mat,1)*3*length(measure_flags),1]);
pm_vec_FDR = mafdr(pm_vec,'BH',1);
pm_FDR = reshape(pm_vec_FDR,[size(VV_mat,1),3,length(measure_flags)]);

% sigstar plots *,**,*** according to p < 0.05, 0.01, 0.001
% modify p-vals to convert *,** to p < 0.05, 0.005
% for labeling in legend to be consistent with other figures
pm = pm_FDR;
pm(pm_FDR>=0.005 & pm_FDR<0.05) = 0.04;
pm(pm_FDR<0.005) = 0.004;

for imeasure = 1:length(measure_flags)
    measure_flag = measure_flags(imeasure);
    V_mat = VV_mat(:,:,imeasure);
    pm0 = pm(:,:,imeasure);

    if measure_flag == 1
        measure_tit = 'Global Graph Measure';
        measure_tag = 'GGM';
    elseif measure_flag == 2
        measure_tit = 'Local Graph Measure';
        measure_tag = 'LGM';
    end

    figure_name = ['Vector Norms of ',measure_tit,' Distance Features'];

    fig = figure('Position',[20,50,1500,350]);
    for isel_feat = 1:size(V_mat,1)
        mz_Dmat0 = V_mat{isel_feat,1};
        sib_Dmat0 = V_mat{isel_feat,2};
        ur_Dmat0 = V_mat{isel_feat,3};
        cat_vec = [ones([1,length(mz_Dmat0)]),2*ones([1,length(sib_Dmat0)]),3*ones([1,length(ur_Dmat0)])];
        data_vec = [mz_Dmat0, sib_Dmat0, ur_Dmat0];

        pFDR_vec = pm0(isel_feat,:);

        subplot(1,size(V_mat,1),isel_feat)
        boxchart(cat_vec,data_vec,'MarkerStyle','.')
        xlim([0.5,3.5])
        xticks(1:length(group_labels))
        set(gca,'XTickLabel',group_labels,'XTickLabelRotation',0)
        title(fig_tits{isel_feat})
        set(gca,'FontSize',12)
        ylims = get(gca,'YLim');

        H = sigstar(groups,pFDR_vec);
        ylim([ylims(1),max(data_vec)+(max(data_vec)-ylims(1))/4.5])
        if strcmp(measure_tag,'GGM')
            ylims2 = get(gca,'YLim');
            yticks(ylims2(1):2*round(diff(ylims2)/12):ylims2(2))
        end
        for icomp = 1:length(groups)
            set(H(icomp,2),'FontSize',14)
        end

    end
    sgtitle(figure_name)

    dpi = '1200';
    img_outfile = [outdir,'\',figure_name];
    img_outfile2 = [outdir,'\',figure_name,'_',dpi,'dpi'];
    save_figure_hp(fig,save_flag,img_outfile2,dpi,0,img_outfile)
end
