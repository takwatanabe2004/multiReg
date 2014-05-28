posx = round(ROI.nx/2)
posy = round(ROI.ny/2)
posz = round(ROI.nz/2)

slicexg = im2uint8(squeeze(refRoi(posx,:,:)))';
sliceyg = im2uint8(squeeze(refRoi(:,posy,:)))';
slicezg = im2uint8(squeeze(refRoi(:,:,posz)))';
slicex = ind2rgb(slicexg,gray(256));
slicey = ind2rgb(sliceyg,gray(256));
slicez = ind2rgb(slicezg,gray(256));

slicex_x = [posx posx; ...
            posx posx]
slicex_y = [0 (ROI.ny-1); ...
            0 (ROI.ny-1)]; 
slicex_z = [0 0; ...
           (ROI.nz-1) (ROI.nz-1)];

slicey_x = [0 (ROI.nx-1); ...
            0 (ROI.nx-1)]; 
slicey_y = [posy posy; ...
            posy posy]; 
slicey_z = [0 0; ...
           (ROI.nz-1) (ROI.nz-1)];

slicez_x = [0 (ROI.nx-1); ...
            0 (ROI.nx-1)]; 
slicez_y = [0 0; ...
            (ROI.ny-1) (ROI.ny-1)]; 
slicez_z = [posz posz; ...
            posz posz];

figure
surface(slicex_x,slicex_y,slicex_z, slicex,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1); hold on
surface(slicey_x,slicey_y,slicey_z, slicey,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1);
surface(slicez_x,slicez_y,slicez_z, slicez,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct','FaceAlpha',1);
rotate3d on
