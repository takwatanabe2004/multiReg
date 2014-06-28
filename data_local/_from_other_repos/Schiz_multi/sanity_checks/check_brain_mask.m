%% save_brain_mask 
% (05/06/2014)
%=========================================================================%
% - sanity check on the brain mask generated from save_brain_mask.m
%=========================================================================%
%%
clear
purge
dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

mask=load_nii('brain_mask.nii')
view_nii(mask)

Subj1_gm = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii'));
Subj1_wm = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp2mprage.nii'));
Subj1_cf = load_nii(strcat(dir_head,'0040000',dir_tail,'m0wrp3mprage.nii'));
Subj1_mask=Subj1_gm;
Subj1_mask.img = (Subj1_gm.img~=0 | Subj1_wm.img~=0 | Subj1_cf.img~=0);

Subj2_gm = load_nii(strcat(dir_head,'0040001',dir_tail,'m0wrp1mprage.nii'));
Subj2_wm = load_nii(strcat(dir_head,'0040001',dir_tail,'m0wrp2mprage.nii'));
Subj2_cf = load_nii(strcat(dir_head,'0040001',dir_tail,'m0wrp3mprage.nii'));
Subj2_mask=Subj2_gm;
Subj2_mask.img = (Subj2_gm.img~=0 | Subj2_wm.img~=0 | Subj2_cf.img~=0);

view_nii(Subj1_mask)
view_nii(Subj2_mask)

view_nii(Subj1_gm)
view_nii(Subj2_gm)