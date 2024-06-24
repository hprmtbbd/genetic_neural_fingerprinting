function [ix_list,ix_title] = atlas_Brainnetome_findLabelIndices(tissuelabel)
func = @(label,labels) cell2mat(cellfun(@(x) strfind(label,x),labels,'UniformOutput',0));

ix_allstruc = 1:length(tissuelabel);
ix_cortical = []; 
ix_subcortical = []; 

ix_allstruc_r = []; ix_allstruc_l = [];
ix_frontal = []; ix_temporal = []; ix_parietal = [];
ix_insular = []; ix_limbic = []; ix_occipital = [];
ix_amyg = []; ix_hipp = []; ix_basal = []; ix_thalamus = [];
for i = 1:length(tissuelabel)
    label = tissuelabel{i};
    if strfind(label,'Right')
        ix_allstruc_r = [ix_allstruc_r,i];
    elseif strfind(label,'Left')
        ix_allstruc_l = [ix_allstruc_l,i];
    end
    
    if any(func(label,{'SFG','MFG','IFG','OrG','PrG','PCL'}))
        ix_frontal = [ix_frontal,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'STG','MTG','ITG','FuG','PhG','pSTS'}))
        ix_temporal = [ix_temporal,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'SPL','IPL','Pcun','PoG'}))
        ix_parietal = [ix_parietal,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'INS'}))
        ix_insular = [ix_insular,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'CG'}))
        ix_limbic = [ix_limbic,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'MVOcC','LOcC'}))
        ix_occipital = [ix_occipital,i];
        ix_cortical = [ix_cortical,i];
    elseif any(func(label,{'Amyg'}))
        ix_amyg = [ix_amyg,i];
        ix_subcortical = [ix_subcortical,i];
    elseif any(func(label,{'Hipp'}))
        ix_hipp = [ix_hipp,i];
        ix_subcortical = [ix_subcortical,i];
    elseif any(func(label,{'BG'}))
        ix_basal = [ix_basal,i];
        ix_subcortical = [ix_subcortical,i];
    elseif any(func(label,{'Tha'}))
        ix_thalamus = [ix_thalamus,i];
        ix_subcortical = [ix_subcortical,i];
    end
end
names = {'Cortical'; 'Subcortical'; 'Frontal'; 'Temporal'; 'Parietal'; 'Insular'; 'Limbic'; 'Occipital';...
    'Amygdala'; 'Hippocampus'; 'Basal Ganglia'; 'Thalamus'};
ix_list_all = {ix_cortical; ix_subcortical; ix_frontal; ix_temporal; ix_parietal; ix_insular; ix_limbic; ix_occipital;...
    ix_amyg; ix_hipp; ix_basal; ix_thalamus};

ix_list_r = {}; ix_list_l = {};
for i = 1:length(ix_list_all)
    ix_temp = ix_list_all{i};
    ix_list_r{i,1} = ix_allstruc_r(ismember(ix_allstruc_r,ix_temp));
    ix_list_l{i,1} = ix_allstruc_l(ismember(ix_allstruc_l,ix_temp));
end

ix_title = {'Name','All_Index','Right_Index','Left_Index'};
ix_list = [{'AllStruc', ix_allstruc, ix_allstruc_r, ix_allstruc_l}; names,ix_list_all,ix_list_r,ix_list_l];
end

