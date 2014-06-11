clear
purge

fsave=true;

parcellation = 'WashU'; % {'Grid326','Grid1068','WashU'}

outPath=[fileparts(mfilename('fullpath')),'/designMatrix_FC_',parcellation,'.mat'];
outVars={'X','SubjList','age','sex','mFileName','timeStamp'};
% return
%%
path_head = '/net/pizza/NKI-RS/FirstLevel/0'; %<- 0 padded to subject ID to match file-directory convention on freewill (6/10/2014)
path_tail = ['/Concat_645_1400/s6_',parcellation,'_p20f0b_p35mask/s6_',parcellation,'_p20f0b_p35mask_roiTC.mat'];

load('/home/slab/users/takanori/multiReg/data_local/csv/MDF_main_info.mat',...
    'SubjList','age','sex')
switch parcellation
    case 'Grid326'
        numNodes=326;        
    case 'Grid1068'
        numNodes=1068;
    case 'WashU'
        numNodes=264;
end
p=nchoosek(numNodes,2);
n = length(SubjList);
%% save design matrix
X = zeros(n,p);
for idx=1:n
    idx
    load(strcat(path_head,SubjList{idx},path_tail),'roiTC');
%     return
%     keyboard
    X(idx,:) =tak_dvec(corr(roiTC))';
    clear roiTC
end
%%
timeStamp=tak_timestamp;
mFileName=mfilename;

if fsave
    save(outPath,outVars{:})
end