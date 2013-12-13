function [plane] = FindTheLid(Image)

% 1. Read the image (NIFTI format)
[VentrNii] = load_nii(Image);
    % Get the matrix of the 3D image:
    Ventr = VentrNii.img;
    % Get the information about image spacing and origin:
    origin = [VentrNii.hdr.hist.qoffset_x VentrNii.hdr.hist.qoffset_y VentrNii.hdr.hist.qoffset_z];
    hPH.origin   = reshape(origin,1,3);
    hPH.spacing  = reshape(VentrNii.hdr.dime.pixdim(2:4),1,3);  
    hPH.dim      = reshape(VentrNii.hdr.dime.dim(2:4),1,3);

% build an insosurface and visualise it:
show_segment_surface(Ventr,hPH.origin,hPH.spacing);

% 2. Set of operations to find the lid
% 2.1: find the blood pools (RV and LV):
    
[R,S,dimensions,LVpool,RVpool] = cav_initial_alignment(Ventr,hPH,'BiV',[]);

% SCRIPT WORKING UP TO HERE!

% 2.2: 
Lid = LVdilated;
se = ones(3,3,3);
LVdilated = imdilate(LVpool,se);
ClosedVentricles = RVpool & LVpool & Ventr; 

% 2.3: find the plane that fits the 3D points you have found, use a
% combination of:
points = getPoints_fromBinary
FitFlatPlane2Points(points);



