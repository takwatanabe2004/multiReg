%%

% we want to grab data from slice Z=22

p = load('/net/data4/Schiz_COBRE/FirstLevel/0040008/347rois_Censor_FDpoint5_0f0bExclude/347rois_Censor_FDpoint5_0f0bExclude_parameters.mat');

coord = p.parameters.rois.mni.coordinates;

nodemask = coord(:,3) == 18;

edgemask = zeros(347);

coord_used = coord(nodemask,:);

for i=1:347
    for j=i:347
        if nodemask(i) & nodemask(j)
            edgemask(i,j)=1;
        end
    end
end

edgemaskflat = logical(mc_flatten_upper_triangle(edgemask));

%% loop over data and extract edgemasked values
SubjDir = {   
   '0040000'	, +1	;
%    '0040001'	, +1	FALSE
   '0040002'	, +1	;
   '0040003'	, +1	;
   '0040004'	, +1	;
%    '0040005'	, +1	FALSE
   '0040006'	, +1	;
%    '0040007'	, +1	FALSE
   '0040008'	, +1	;
%    '0040009'	, +1	FALSE
   '0040010'	, +1	;
%    '0040011'	, +1	FALSE
   '0040012'	, +1	;
   '0040013'	, -1	;
   '0040014'	, -1	;
   '0040015'	, +1	;
%    '0040016'	, +1	FALSE
   '0040017'	, -1	;
   '0040018'	, -1	;
   '0040019'	, -1	;
   '0040020'	, -1	;
   '0040021'	, +1	;
   '0040022'	, +1	;
   '0040023'	, -1	;
%    '0040024'	, -1	FALSE
   '0040025'	, +1	;
%    '0040026'	, -1	FALSE
   '0040027'	, -1	;
   '0040028'	, +1	;
   '0040029'	, +1	;
%    '0040030'	, -1	FALSE
%    '0040031'	, -1	FALSE
   '0040032'	, +1	;
   '0040033'	, -1	;
   '0040034'	, +1	;
   '0040035'	, -1	;
%    '0040036'	, -1	FALSE
   '0040037'	, +1	;
   '0040038'	, -1	;
%    '0040039'	, +1	FALSE
%    '0040040'	, +1	FALSE
   '0040041'	, +1	;
   '0040042'	, +1	;
   '0040043'	, -1	;
   '0040044'	, +1	;
   '0040045'	, -1	;
   '0040046'	, +1	;
   '0040047'	, +1	;
   '0040048'	, -1	;
   '0040049'	, +1	;
   '0040050'	, -1	;
   '0040051'	, -1	;
   '0040052'	, -1	;
   '0040053'	, -1	;
   '0040054'	, -1	;
   '0040055'	, -1	;
   '0040056'	, -1	;
   '0040057'	, -1	;
   '0040058'	, -1	;
   '0040059'	, +1	;
%    '0040060'	, +1	FALSE
   '0040061'	, -1	;
   '0040062'	, -1	;
   '0040063'	, -1	;
%    '0040064'	, +1	FALSE
   '0040065'	, -1	;
   '0040066'	, -1	;
   '0040067'	, -1	;
   '0040068'	, -1	;
   '0040069'	, -1	;
%    '0040070'	NA	FALSE
%    '0040071'	, +1	FALSE
   '0040072'	, +1	;
%    '0040073'	, +1	FALSE
   '0040074'	, -1	;
%    '0040075'	, +1	FALSE
   '0040076'	, -1	;
   '0040077'	, +1	;
   '0040078'	, +1	;
   '0040079'	, +1	;
   '0040080'	, +1	;
   '0040081'	, +1	;
   '0040082'	, +1	;
%    '0040083'	NA	FALSE
%    '0040084'	, +1	FALSE
   '0040085'	, +1	;
   '0040086'	, -1	;
   '0040087'	, -1	;
%    '0040088'	, +1	FALSE
%    '0040089'	, +1	FALSE
   '0040090'	, -1	;
   '0040091'	, -1	;
   '0040092'	, +1	;
   '0040093'	, -1	;
   '0040094'	, +1	;
   '0040095'	, -1	;
   '0040096'	, +1	;
%    '0040097'	, +1	FALSE
   '0040098'	, +1	;
   '0040099'	, +1	;
   '0040100'	, +1	;
   '0040101'	, +1	;
   '0040102'	, -1	;
   '0040103'	, +1	;
   '0040104'	, -1	;
%    '0040105'	, +1	FALSE
   '0040106'	, +1	;
   '0040107'	, -1	;
   '0040108'	, +1	;
   '0040109'	, +1	;
   '0040110'	, +1	;
%    '0040111'	, -1	FALSE
   '0040112'	, +1	;
   '0040113'	, -1	;
   '0040114'	, -1	;
   '0040115'	, -1	;
   '0040116'	, -1	;
   '0040117'	, +1	;
   '0040118'	, -1	;
   '0040119'	, -1	;
   '0040120'	, -1	;
   '0040121'	, -1	;
   '0040122'	, +1	;
   '0040123'	, -1	;
   '0040124'	, -1	;
   '0040125'	, -1	;
   '0040126'	, +1	;
   '0040127'	, -1	;
   '0040128'	, -1	;
   '0040129'	, -1	;
   '0040130'	, -1	;
   '0040131'	, -1	;
   '0040132'	, +1	;
   '0040133'	, +1	;
   '0040134'	, -1	;
   '0040135'	, -1	;
%    '0040136'	, -1	FALSE
   '0040137'	, +1	;
   '0040138'	, -1	;
   '0040139'	, -1	;
   '0040140'	, -1	;
   '0040141'	, -1	;
   '0040142'	, +1	;
   '0040143'	, +1	;
   '0040144'	, -1	;
   '0040145'	, +1	;
   '0040146'	, -1	;
   '0040147'	, -1	;
};

SubjDirHC = SubjDir((cell2mat(SubjDir(:,2)))==-1,:);

nSub = size(SubjDirHC,1);
nNode = sum(nodemask);

mini = zeros(nNode);
nMiniflat = numel(mc_flatten_upper_triangle(mini));

temp = '/net/data4/Schiz_COBRE/FirstLevel/[Subject]/347rois_Censor_FDpoint5_0f0bExclude/347rois_Censor_FDpoint5_0f0bExclude_corr.mat';


data = zeros(nSub,nMiniflat);

for i = 1:size(SubjDirHC,1)
    Subject = SubjDirHC{i,1};
    path = mc_GenPath(temp);
    conn = load(path);
    rflat = mc_flatten_upper_triangle(conn.rMatrix);
    data(i,:) = rflat(edgemaskflat);
end

data = mc_fisherz(data);

ConnMeanFlat = mean(data,1);
ConnVarFlat = var(data,1);

ConnMean = mc_unflatten_upper_triangle(ConnMeanFlat);
ConnVar = mc_unflatten_upper_triangle(ConnVarFlat);

cd /net/data4/Schiz_COBRE/Scripts/slab/Tak_ADMM/
save('Results.mat','ConnMean','ConnVar','edgemask','coord_used','coord','nodemask')