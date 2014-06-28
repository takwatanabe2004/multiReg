clear
purge

dir_head = 'C:\Users\takanori\Desktop\Subjects/';
dir_tail = '/session_1/rest_1/coReg/';

dir_nii1=strcat(dir_head,'0040000',dir_tail,'m0wrp1mprage.nii')
nii1 = load_nii(dir_nii1);
dir_nii2=strcat(dir_head,'0040001',dir_tail,'m0wrp1mprage.nii')
nii2 = load_nii(dir_nii2);

view_nii(nii1)
view_nii(nii2)
sum(nii1.img(:)==0)
sum(nii2.img(:)==0)
% return
%%
% purge
nii1.img=nii1.img~=0;
view_nii(nii1)
nii2.img=nii2.img~=0;
view_nii(nii2)
% wm_mask = load_nii('csf_mask.nii');
% view_nii(wm_mask)
return
%%
% nii.img = imfill(nii.img,[3 3 3],26);
% view_nii(nii)
%%
img_vec=img(:);
img_vec(img_vec==0)=[];
% hist(img_vec,1000)

perc_zeroes=sum(img(:)==0)/numel(img)*100
perc_isnan=sum(isnan(img(:)))/numel(img)*100