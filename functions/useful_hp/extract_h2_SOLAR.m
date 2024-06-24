function [hval,pval,herror,covar_var,kurt,kurt_score] = extract_h2_SOLAR(indir0)

txt_file = fullfile(indir0,'polygenic.out');
fD = fopen(txt_file);
txt = textscan(fD,'%s');
fclose(fD);

hind = find(strcmp('H2r',txt{1}));
hval = str2double(txt{1}{hind(1)+2});
pval = str2double(txt{1}{hind(1)+5});
if length(hind) < 2
    herror = nan;
else
    herror = str2double(txt{1}{hind(2)+3});
end

varind = find(strcmp('Covariates',txt{1}));
covar_var = str2double(txt{1}{varind(1)+2});

kurtind = find(strcmp('Kurtosis',txt{1}));
kurt = str2double(txt{1}{kurtind(1)+2});
kurt_score = txt{1}{kurtind(1)+4};

end
