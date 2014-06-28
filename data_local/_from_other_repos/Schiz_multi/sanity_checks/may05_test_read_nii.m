clear
purge
nii = load_nii('m0wrp1mprage.nii');
img=nii.img;
view_nii(nii)
nii.img=img~=0;
view_nii(nii)

% wm_mask = load_nii('csf_mask.nii');
% view_nii(wm_mask)
% return
%%
img_vec=img(:);
img_vec(img_vec==0)=[];
% hist(img_vec,1000)

perc_zeroes=sum(img(:)==0)/numel(img)*100
perc_isnan=sum(isnan(img(:)))/numel(img)*100