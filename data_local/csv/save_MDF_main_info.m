%% save_MDF_vars_of_interest
% (06/10/2014)
%=========================================================================%
% save variables of interest from MDF_NKIRS_ICA.csv
% - (age, sex, subjects, inclusion mask)
%=========================================================================%
%%
clear
purge

fsave=true;

outVars={'SubjList','age','sex','timeStamp','mFileName'}
outPath=[fileparts(mfilename('fullpath')),'/MDF_main_info.mat']

getMainInfoMDF_NKIRS

%% get subject inclusion mask
% R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),mask); % Find non-numeric cells

%=========================================================================%
% get subject inclusion mask
%=========================================================================%
maske = cellfun(@(x) isequal('TRUE',x),Include_645_1400_3min_r4_1); 
nSubs = sum(maske);

%=========================================================================%
% get subject id list (format n cell of strings)
%=========================================================================%
SubjList=num2cell(Subject(maske));
for i=1:nSubs
    SubjList{i}=num2str(SubjList{i});
end

%=========================================================================%
% age of included subjects
%=========================================================================%
age = Age(maske);

%=========================================================================% 
% sex: +1 = male, -1 = female
%=========================================================================%
sex = Sex(maske);
sex = cellfun(@(x) isequal('M',x),sex); 
sex=double(sex);
sex(sex==0)=-1;
%% save
mFileName=mfilename;
timeStamp=tak_timestamp;
if fsave
    save(outPath,outVars{:})
end