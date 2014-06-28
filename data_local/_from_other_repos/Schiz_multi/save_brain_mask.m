%% save_brain_mask
% (05/06/2014)
%=========================================================================%
% (05/08/2014)
% - SCRIPT OBSOLETE!!! "OFFICIAL BRAIN MASK" OBTAINED FROM MIKE!
%-------------------------------------------------------------------------%
% - save a mask of the brain (output: brain_mask.nii)
%=========================================================================%
%%
% clear
% purge
% 
% SubjDir = S_subjDir_Schiz_COBRE_censor;
% SubjList = SubjDir(:,1);
% % return
% %%
% dir_head = 'C:\Users\takanori\Desktop\Subjects/';
% dir_tail = '/session_1/rest_1/coReg/';
% 
% %% initialize first nii struct
% master = load_nii(strcat(dir_head,SubjDir{1},dir_tail,'m0wrp1mprage.nii'))
% master.img = false(size(master.img));
% % return
% %% loop through all masks
% for idx=1:length(SubjList)
%     idx
%     Subj_gm = load_nii(strcat(dir_head,SubjDir{idx},dir_tail,'m0wrp1mprage.nii'));
%     Subj_wm = load_nii(strcat(dir_head,SubjDir{idx},dir_tail,'m0wrp2mprage.nii'));
%     Subj_cf = load_nii(strcat(dir_head,SubjDir{idx},dir_tail,'m0wrp3mprage.nii'));
%     
%     %% set brain region value to 1 (0 for non brain region...which has 0 or nan voxel values)
%     mask_gm = (Subj_gm.img~=0) & ~isnan(Subj_gm.img);
%     mask_wm = (Subj_wm.img~=0) & ~isnan(Subj_wm.img);
%     mask_cf = (Subj_cf.img~=0) & ~isnan(Subj_cf.img);
%     
%     master.img = (master.img | mask_gm | mask_wm | mask_cf);
% end
% view_nii(master)
% 
% % return
% %%
% master.hdr.dime.datatype=1; 
% save_nii(master,'brain_mask.nii')
% 
% %%
