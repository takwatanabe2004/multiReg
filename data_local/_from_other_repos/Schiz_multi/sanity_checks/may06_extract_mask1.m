%% may06_extract_mask1
% (05/06/2014)
%-------------------------------------------------------------------------%
% Problem: subject-to-subject variation exist in the mask...
%-------------------------------------------------------------------------%
%%
clear
purge

dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

Subj1 = '0040000';
Subj1_gm=strcat(dir_head,Subj1,dir_tail,'m0wrp1mprage.nii');
Subj1_wm=strcat(dir_head,Subj1,dir_tail,'m0wrp2mprage.nii');
Subj1_cf=strcat(dir_head,Subj1,dir_tail,'m0wrp3mprage.nii');
Subj1_gm = load_nii(Subj1_gm);
Subj1_wm = load_nii(Subj1_wm);
Subj1_cf = load_nii(Subj1_cf);

Subj2 = '0040001';
Subj2_gm=strcat(dir_head,Subj2,dir_tail,'m0wrp1mprage.nii');
Subj2_wm=strcat(dir_head,Subj2,dir_tail,'m0wrp2mprage.nii');
Subj2_cf=strcat(dir_head,Subj2,dir_tail,'m0wrp3mprage.nii');
Subj2_gm = load_nii(Subj2_gm);
Subj2_wm = load_nii(Subj2_wm);
Subj2_cf = load_nii(Subj2_cf);
%%
view_nii(Subj1_gm)
view_nii(Subj2_gm)
% view_nii(nii2)
% view_nii(nii3)
% return
%%
% purge
Subj1_mask=Subj1_gm;
Subj1_mask.img = (Subj1_gm.img~=0 | Subj1_wm.img~=0 | Subj1_cf.img~=0);

Subj2_mask=Subj2_gm;
Subj2_mask.img = (Subj2_gm.img~=0 | Subj2_wm.img~=0 | Subj2_cf.img~=0);

view_nii(Subj1_mask)
view_nii(Subj2_mask)
% return
%%
sum(Subj1_mask.img(:)~=0)
sum(Subj2_mask.img(:)~=0)