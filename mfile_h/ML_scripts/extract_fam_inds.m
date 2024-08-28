function [imz,jmz,isib,jsib,nmz,nsib,fam_locs00,fam_locs] = extract_fam_inds(phi,famid,nsubj,data_flag)

phi0 = phi;
phi0(1:nsubj+1:end) = 0;

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

end