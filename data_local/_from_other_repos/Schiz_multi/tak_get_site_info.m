function [SubjDir,ylabel,nSubs,nSubsDS,nSubsHC] = tak_get_site_info(options)
% [SubjDir,ylabel,nSubs,nSubsDS,nSubsHC] = tak_get_site_info(options)
%|------------------------------------------------------------------------------|%
%| Returns various information about a specified 'site' from a specific 'disorder'
%| study.
%| 
%| The workspace may get cluttered using this code, so consider using 
%| tak_get_site_info_struct.m, which does the exact same thing, but returns 
%| the output as a 'struct'
%|------------------------------------------------------------------------------|%
%| INPUT
%| options: struct containing the following fields
%| -    disorder: disorder-type - {'Autism', 'ADHD', 'Schiz_COBRE'} 
%| -        site: select site
%| -      motion: {'censor','nocensor'}
%| 
%|  The following 'site' is valid (depends on options.disorder)
%|  ADHD   = {'KKI', 'NeuroIMAGE', 'NYU', 'OHSU', 'Pecking'}
%|  Autism = {'KKI', 'Leuven_1','Leuven_2', 'NYU', 'OHSU', 'Olin', 'Pitt', ...
%|            'SBL','SDSU','Trinity','UCLA_1','UCLA_2','UM_1','UM_2','USM','Yale'}; 
%|  Schiz_COBRE = {''} (empty bracket)
%|------------------------------------------------------------------------------|%
%| OUTPUT
%| - SubjDir: (nSubs x 2) cell array {col1=subjID,col2=class label}
%| -  ylabel: (n x 1) class label vector with entries +/-1
%| -   nSubs: numbers of subjects from the 'site'
%| - nSubsDS: numbers of subjects diagnosed as 'disorder group'
%| - nSubsHC: numbers of subjects diagnosed as 'healty control'
%|------------------------------------------------------------------------------|%
%| 03/05/2013 -> added support for dataset for three disorders.
%| 03/25/2013 -> input is now a struct with three fields
%|               (1) disorder (2) site (3) motion sensor option
%|               Currently, Autism is disabled.
%| 05/18/2013 -> Autism enabled
%%
disorder = options.disorder;
site     = options.site;
motion   = options.motion;

switch disorder
    case 'ADHD'
        SubjDir = tak_subjDir_ADHD_FirstLevel(site,motion);
    case 'Autism'
        SubjDir = tak_subjDir_Autism_FirstLevel(site,motion);
    case 'Schiz_COBRE'
        SubjDir = tak_subjDir_Schiz_COBRE_FirstLevel(motion);
end
%%
nSubs = length(SubjDir);
nSubsDS = sum(cell2mat(SubjDir(:,2))==+1);
nSubsHC = sum(cell2mat(SubjDir(:,2))==-1);

ylabel = cell2mat(SubjDir(:,2));  