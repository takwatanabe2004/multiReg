clear
purge
fid=fopen('MDF_NKIRS_ICA.csv');
out = textscan(fid,'%s%f$f','delimiter',',');
fclose(fid);

date = datevec(out{1});
col1 = out{2};
col2 = out{3};