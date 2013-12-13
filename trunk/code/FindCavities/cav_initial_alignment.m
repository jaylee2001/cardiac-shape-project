% Initial alignment using ventricle cavities.
% Written by Matthew Sinclair
% Created: 21/12/11
% ------------------------------------------------------------------------
%   DESCRIPTION:
%   This script obtains a coordinate system based on the LV (and RV)
%   cavit(y/ies) of a segmentation for a more robust initial alignment. The
%   first principal component of the LV cavity is used as the z-axis (Vz).
%   The component of the vector between the COM of the LV and RV cavities
%   perpendicular to Vz is the y-axis (Vy).
%
%   NOTE: the opensource toolbox iso2mesh is used for the function s2v (at
%   line 80), download from:
%   http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Download
% ------------------------------------------------------------------------
%   INPUTS:
%   1 - Name of segmentation ('--.vtk')
%   2 - Type of segmentation ('LV' or 'BiV')
%   3 - optional debug flag (any integer>0) for visualising intermediate steps
% ------------------------------------------------------------------------
%   OUTPUTS:
%   1 - Rotation vector - R(:,1) is the main inertial axis
%   2 - Scale of LV blood pool
% ------------------------------------------------------------------------
%   TO DO:
%   - Implement feature extraction for tailoring template
% ------------------------------------------------------------------------
%   EXAMPLE USAGE:
%   [R] = cav_initial_alignment('case15VentriclesCleanSmall.vtk','BiV')
%   for debug (figures), add 3rd argument: 1
% ------------------------------------------------------------------------
%
% Revision control:
% 04/09/2013: removal of bOffsetXangleinBiValignment when horizontal
% 15/05/2013: add options to check the correctness of the right to left
% direction (second axis) and upside-down (main axis) based on the
% BasalPlane (options: RVdirection and LVbasalPoint/RVbasalPoint)
% 17/12/2012: add the default bSparseData = 0 at init of funciton
% 01/10/2012: possibility to prescribe the cardiac orientation (R)
% 17/04/2012: possibility to extract cavity dimensions
% 02/04/12: addition of option "bTemplate" in order to speed up the process
% and add the robustness to the definition of the orientation of the LV,
% when there is revolutionary symetry
% 01/03/12: Correction of the scale of whole segmentation, given the fact
% that the orientation of the shape is not any more the result of its PCA,
% whearas the orientation calculated of the LV.
% 27/02/12: reverted calculation of the scale to whole segmentation (M.Sinclair)
% 18/01/12: introduction of an optional flag indicating that the binary
% mask has its basal plane horizontal (P.Lamata)
% 26/12/11: change of the input of the function (P.Lamata)
% 21/12/11: calculation of the scale of the LV blood pool (P.Lamata)

function [R,S,dimensions,Lpool,Rpool] = cav_initial_alignment(seg_im,seg_hd,seg_type,options)

%==========================================================================
% Default options:

% In BiV shapes, the LV main axis is usually tilted with respect to the
% base plane, and this is heuristically fixed to:
OffsetXangle = 10; % in degrees
switch seg_type
    case 'LV'
        bOffsetXangleinBiValignment = 0;
    otherwise
        bOffsetXangleinBiValignment = 1;
end
S = ones(1,3);
debug = 0;
bDebugRotationMatrix = 0;
bHorizontalBase = 0;

bLVdirGiven     = 0;
bRVdirectionGiven = 0;
bRVdirGivenHorizontalPlane     = 0;
bRVpoolLabel    = 0;

bRotationGiven  = 0;
bReferenceRotationGiven = 0;
bTemplate       = 0;
bDimensions     = 0;
bRScalculated   = 0;
bSparseData     = 0;
seg_name = '';
SegmentationLabel = [];
Lpool = [];
Rpool = [];
%==========================================================================

if nargin==4
    if isfield(options,'OffsetXangle'), OffsetXangle = options.OffsetXangle; end
    if isfield(options,'debug'), debug = options.debug; end
    if isfield(options,'bHorizontalBase'), bHorizontalBase = options.bHorizontalBase; end
    if isfield(options,'bTemplate'), bTemplate = options.bTemplate; end
    if isfield(options,'bDimensions'), bDimensions = options.bDimensions; end
    if isfield(options,'RotationMatrix'), R = options.RotationMatrix; bRotationGiven = 1; end
    if isfield(options,'SegmentationLabel'), SegmentationLabel = options.SegmentationLabel; end
    if isfield(options,'LVdir'), 
        % directon of the main axis of inertia:
        LVdir = options.LVdir;
        bLVdirGiven = 1;
    end
    if isfield(options,'RVdir'),
        % This is the direction in an horizontal plane (as in a short axis
        % stack of images)
        RVdir = options.RVdir;         
        if numel(RVdir)==2 
            bRVdirGivenHorizontalPlane = 1;
        else
            fprintf('ERROR! The parameter RVdir is invalid ignored in cav_initial_alignment\n');
            RVdir
        end
    end
    if isfield(options,'RVdirection'),
        % This is the direction in 3D
        RVdirection = options.RVdirection;         
        optionsGetRotMatrix.RVdirection = RVdirection;
        if numel(RVdirection)==3 
            bRVdirectionGiven = 1;
        else
            fprintf('ERROR! The parameter RVdirection is invalid ignored in cav_initial_alignment\n');
            RVdirection
        end
    end     
    if isfield(options,'RVpoolLabel'),  
        bRVpoolLabel = 1;
        RVpoolLabel  = options.RVpoolLabel;  
        RVbp = zeros(size(seg_im));
        RVbp(seg_im == RVpoolLabel) = 1;
        optionsGetRotMatrix.RVbp = RVbp;
    end
    if isfield(options,'LVbasalPoint'), optionsGetRotMatrix.LVbasalPoint = options.LVbasalPoint; end
    if isfield(options,'RVbasalPoint'), optionsGetRotMatrix.RVbasalPoint = options.RVbasalPoint; end
end

if(bHorizontalBase)
    bOffsetXangleinBiValignment = 0;
end
    

% Read input file
% Change by P.Lamata, 21/12/2011:
%[seg_im seg_hd] = io_ReadMedicalImage(seg_name);
% if strcmp(seg_name(end-3:end),'.vtk')
%     read_image_vtk2(seg_name);
% elseif strcmp(seg_name(end-4:end),'.gipl')
%     [seg_im seg_hd] = ReadData3D(seg_name);
% end

if isfield(options,'SegmentationLabel')
	seg_im = GetSegmentationByLabels(seg_im,SegmentationLabel);
%     seg_im = boolean(seg_im);
    seg_im = logical(seg_im);   %Trying to convert to boolean matrix? G.Gonzalez 24/09/13
else
    % It might be that the image has more than 2 labels, it is not a
    % binary:
    nLabels = unique(seg_im);
    if numel(nLabels)>2     % label=0, blank space; label>0, segmentation
        % It might be that the image is a grayscale. Decision is based on
        % the number of labels:
        if numel(nLabels) > 50
            % This is a grayscale image. Take all values above 0 as
            % myocardium:
            seg_im(seg_im>0) = 1;
        else
            % Find the biggest region segmented as the myocardium:
            seg_im = find_myo(seg_im,seg_hd,seg_name,debug);
        end
    end
end

% TODO: Find suitable truncation (May 2013: done before calling this
% function, see heartgen2.m)
% P.Lamata: 13/06/2012 remove unconnected bits to prevent segmentation
% noise crash the detection of blood pools
seg_im = RemoveSmallestComponents(seg_im,10000,26,1);

% convert to points
seg_p = getPoints_fromBinary(seg_im,seg_hd.origin,seg_hd.spacing);

% Check if rotation is already given by RV and LV directions:
if(bLVdirGiven) && (bRVdirectionGiven)
    fprintf('Rotation defined by LV and RV directions (input in parameters)\n')
    % LV is V1, and RV is V2, so we can assemble the rotation matrix:
    v1 = LVdir;
    v2 = RVdirection;
    PCAt = PCAtransformation;
    Rref = PCAt.BuildRotMatrixFromTwoVectors(v1,v2);    
    bReferenceRotationGiven = 1;
end


% The rotation matrix could be given, and only needing the scale:
if(bRotationGiven)
    fprintf('Rotation of shape given as input, only need to compute scale\n');
    S = GetScaleFromPointsAndAxis(seg_p,R);
    bRScalculated = 1;
else    
    if(bTemplate)
        paramsTemplateAlignment.bOffsetXangleinBiValignment = bOffsetXangleinBiValignment;
        paramsTemplateAlignment.OffsetXangle = OffsetXangle;
        fprintf('\n ** Computation of template rotation and scale ** \n');
        [R,S] = Template_Initial_Alignment(seg_p,paramsTemplateAlignment); 
        bRScalculated = 1;
    else
        fprintf('\n ** Computation of shape rotation and scale ** \n');
        % Shortcut in case of sparse dataset:
%         if isfield(seg_hd,'OriginalSpacing')
%             % Images loaded from file would have this, if they've gone through
%             % ParseHeader
%             bSparseData = DetectSparseData(seg_hd.OriginalSpacing);
%         else
%             % Images might come from a model, or elsewhere:
%             bSparseData = DetectSparseData(seg_hd.spacing);
%         end       

        % The search for the blood pool was required to identify the right scale of
        % the blood pool!
        % This bit of code is taken back (P.Lamata), since sparse cases otherwise
        % fail in the search of the blood pools, and this option is much quicker!
        if(bSparseData)||bHorizontalBase
            fprintf('Cavity Initial Alignment: base is horizontal (bSparse=%i) (bHorizontalBase=%i)\n',bSparseData,bHorizontalBase)
            % In case of sparse data, the main vector is pointing in the slice
            % direction. In case of a LV, nothing else is needed.
            if strcmp(seg_type,'LV')           
                if(bRVdirGivenHorizontalPlane)
                    [R,S] = Assemble_Horizontal_Initial_Alignment(seg_p,RVdir);
                else
                    [R,S] = Horizontal_Inital_Alignment(seg_p);
                end
                % There is the additional option of the introduction of the
                % RV blood pool label that will rearrange this initial
                % vertical alignment:
                if(~bRVpoolLabel)
                    bRScalculated = 1;
                end
            end
        end
    end
end
if(~bRScalculated)||(bDimensions)
    % convex hull
    seg_ch_k = convhulln(seg_p);

    % Visualise for debug -----------------------------------------------------
    if debug
        figure; show_segment_surface(seg_im,seg_hd.origin,seg_hd.spacing,0.1,0.5);
        hold on; trisurf(seg_ch_k,seg_p(:,1),seg_p(:,2),seg_p(:,3));
        if isfield(seg_hd,'File')
            seg_name = seg_hd.File;
        else seg_name = 'Template';
        end
        title([seg_name ' with convex hull']);
        axis equal;
        input('Press enter to close figure and continue in debug mode.');
        close;
    end
    % -------------------------------------------------------------------------


    % convert convex hull points to shell volume
    dx = min(seg_hd.spacing);   % smallest spacing dimension
    p0=min(seg_p);              % minimum coordinates
    p1=max(seg_p);              % maximum coordinates
    lim = p1-p0;                % difference
    div = min(lim)/dx;          % number of divisions
    seg_ch_im = s2v(seg_p,seg_ch_k,div);    % convert surface points to binary shell


    % convert shell to solid volume
    seg_ch_im_filled = logical(imfill(seg_ch_im,'holes'));

    % if 'holes' fails.
    if numel(find(seg_ch_im_filled)) == numel(find(seg_ch_im))
        bg = logical((imfill(seg_ch_im,[1 1 1])));
        seg_ch_im_filled = (bg==0) + seg_ch_im;
    end

    % Re-bin segmentation binary to match convex hull binary (mask)
    seg_im_rebin = (zeros(size(seg_ch_im_filled,1),size(seg_ch_im_filled,2),size(seg_ch_im_filled,3)));
    seg_im_rebin_temp = getBinary_fromPoints(seg_p,ceil([lim(1)/dx lim(2)/dx lim(3)/dx]),p0,[dx dx dx]);
    seg_im_rebin(1:size(seg_im_rebin_temp,1),1:size(seg_im_rebin_temp,2),1:size(seg_im_rebin_temp,3)) = logical(seg_im_rebin_temp);


    % Visualise for debug -----------------------------------------------------
    if (debug)
        hold on; show_segment_surface(seg_im,seg_hd.origin,seg_hd.spacing,0.1,0.5);
        hold on; show_segment_surface(seg_ch_im_filled,p0,[dx dx dx],0.9,0.5);
        title('Agreement between (rebinned) convex hull (red) and original segmentation (blue)');
        axis equal;
        input('Press enter to close figure and continue in debug mode.');
        close;

        imagesc(squeeze(max(seg_ch_im_filled(:,:,:),[],3)));
        title('Z-MIP of convex hull');
        input('Press enter to compare Z-MIP of convex hull with rebinned segmentation.');
        hold on; imagesc(squeeze(max(seg_im_rebin(:,:,:),[],3)));
        title('Z-MIP of rebinned segmentation');
        input('Press enter to close figure and continue in debug mode.');
        close;
    end
    % % -------------------------------------------------------------------------


    % Obtain Rotation matrix from blood pools
    if bSparseData||bHorizontalBase
        bZMainAxis = 1;
    else
        bZMainAxis = 0;
    end
    optionsGetRotMatrix.bZMainAxis = bZMainAxis;
    if(bLVdirGiven)
        optionsGetRotMatrix.LVdir = LVdir;
    end
    if strcmp(seg_type,'LV')
        LV = find_cavities(seg_ch_im_filled,seg_im_rebin,seg_type,p0,dx,debug);
    elseif strcmp(seg_type,'BiV')
        [LV RV] = find_cavities(seg_ch_im_filled,seg_im_rebin,seg_type,p0,dx,debug);        
        optionsGetRotMatrix.RV = RV;
        Lpool = LV;
        Rpool = RV;
    end
    nVoxelsLV = numel(LV(LV>0));
    if nVoxelsLV<=3
        fprintf('ERROR! Blood pool found is too small (it has %i voxels)!!\n',nVoxelsLV)
        return;
    else
        [Rtemp,S,bSwapLV_RV] = get_rotation_matrix(seg_type, LV, p0, dx, seg_p, optionsGetRotMatrix);
        if(bSwapLV_RV)
            fprintf(' SWAP LV AND RV!\n');
            fprintf('   It looks like the RV had more voxels than the LV\n');
            fprintf('   and that the reference direction helped to correct it\n');
            temp = LV;
            LV = RV;
            RV = temp;
        end
        if(~bRScalculated)
            R = Rtemp;
        end
        if(bReferenceRotationGiven)
            % Compare rotation matrixes
            RotDif = R - Rref;
            epsilon = 1;
            if (sum(abs(RotDif(:)))>epsilon)
                fprintf('Apparent orientation error!\n')
                fprintf(' Reference orientation (from LVdir and RVdir) is:\n')
                fprintf(' %1.2f %1.2f %1.2f\n %1.2f %1.2f %1.2f\n %1.2f %1.2f %1.2f\n',Rref);
                fprintf(' Computed orientation is:\n')
                fprintf(' %1.2f %1.2f %1.2f\n %1.2f %1.2f %1.2f\n %1.2f %1.2f %1.2f\n',R);
            end
        end
        % Visualise for debug -----------------------------------------------------
        if debug
            % transform segmentation to find dimensions for scaling
            %R2 = R;
            %R(:,1) = -R2(:,3);R(:,3) = R2(:,1);
            seg_p_t = (R'*seg_p')';         % transform segmentation
            p_b0=min(seg_p_t);              % minimum coordinates
            p_b1=max(seg_p_t);              % maximum coordinates
            lim_b = p_b1-p_b0;
            seg_im_t = getBinary_fromPoints(seg_p_t,[lim_b(1) lim_b(2) lim_b(3)],p_b0,[1 1 1]);

            LV_p = getPoints_fromBinary(LV,p0,[dx dx dx]);
            LV_p_t = (R'*LV_p')';
            LV_im_t = getBinary_fromPoints(LV_p_t,[lim_b(1) lim_b(2) lim_b(3)],p_b0,[1 1 1]);

            show_segment_surface(seg_im_t,p_b0,[dx dx dx],0.8,0.1);
            show_segment_surface(LV_im_t,p_b0,[dx dx dx],0.3,0.8);
            axis equal;
            view(90,0);
            input('Press enter to close figure and continue in debug mode.');
            close;
        end
    end
    % -------------------------------------------------------------------------
end


Centre = mean(seg_p,1);
if(bDebugRotationMatrix)
    figure('color',[1 1 1]); hold on;
    % Visualise the geometry with the main orientation vectors:
    show_segment_surface(seg_im,seg_hd.origin,seg_hd.spacing);
    color = 'rgb';
    for iC = 1:3
        Vector = R(:,iC);
        P = zeros(3,2);
        for jC=1:3
            P(jC,:) = [Centre(jC) Centre(jC)+Vector(jC)*sqrt(S(jC,jC))];
        end
        plot3(P(1,:),P(2,:),P(3,:),color(iC),'LineWidth',3);
    end
    title('Main axis of shape. 1:red, 2:green, 3:blue');
    axis equal; 
    pause();
end 

% TODO: Obtain Rotation matrix from flat base

dimensions.Centre       = Centre;
dimensions.RVapex       = NaN;
dimensions.RVsizeM      = NaN;
dimensions.RVsizeminor  = NaN;
dimensions.MainLVAxis   = NaN;
dimensions.Length       = NaN;
dimensions.SeptumAngle  = NaN;
dimensions.bBiVw        = NaN;
dimensions.hIm          = seg_hd;
dimensions.WallT        = NaN;
if(bDimensions)    
    if ~exist('LV','var')
        % Need to get the points that are part of the LV and RV for further
        % analysis!
        if strcmp(seg_type,'LV')
            LV = seg_im;
            p0 = seg_hd.origin;
            SamplingResolution = seg_hd.spacing;
        else
            fprintf('Code not prepared yet to analyse dimensions in the conditions used!\n')
        end
    else
        SamplingResolution = [dx dx dx];
    end
    FirstAxis = R(:,1);
    SecondAxis = R(:,2); 
    ThirdAxis = R(:,3);               
    % Wall thickness:
    
    % COMPUTE WALL THICKNESS:
    % Check if this is a sparse case:
    MeanVoxelLength = mean(dimensions.hIm.spacing(1:2));
    MeanVoxelSize = mean(dimensions.hIm.spacing(1:3));
    if DetectSparseData(dimensions.hIm.spacing)
        [WallThickness] = EstimateWallThicknessFromSA(seg_im,seg_type,R);        
        % conversion to mm:
        WallThickness = WallThickness * MeanVoxelLength;
        LengthVoxelSize = dimensions.hIm.spacing(3);
    else
        LengthVoxelSize = MeanVoxelSize;
        % TODO: with a full resolution anatomy
        WallThickness = NaN;
        %[WallThickness] = EstimateWallThicknessFromSA(seg_im,seg_type,R); 
        
    end
    dimensions.WallT = WallThickness;
    
    % Get the reference of the LV:
    LV_p = getPoints_fromBinary(LV,p0,SamplingResolution);       
    LVlength = GetStandardDeviationInDirection(LV_p,FirstAxis); 
    LVminorAxis(1) = GetStandardDeviationInDirection(LV_p,SecondAxis);
    LVminorAxis(2) = GetStandardDeviationInDirection(LV_p,ThirdAxis);
    LVaxis = mean(LVminorAxis);
    
    % The estimation of the length is taking two standard deviations: the
    % scalars here are chosen empirically (shapes were coming too elongated
    % otherwise):
    % Some reference values from Wikipedia:
    % - diameter: 48mm = MainLVAxis
    NumberStdsLength = 4; % before 1.5
    if(bSparseData)&&strcmp(seg_type,'LV')
        % Fine tuned with the cohort of 40 preclampsia cases
        NumberStdsRadius = 5;
    else
        % Fine tuned with the two average anatomies. 
        NumberStdsRadius = 7.5;
    end
    % 12/03/13: removal of the "MeanVoxelSize" factor, units are already in
    % mm!
    dimensions.MainLVAxis  = NumberStdsRadius * LVaxis;% * MeanVoxelSize/dx;
    dimensions.Length = NumberStdsLength * LVlength;% * LengthVoxelSize/dx;
    dimensions.bSparseData = bSparseData;
        
        
    if strcmp(seg_type,'BiV')
        % Septum angle:

        % RV major axis (in the left to right direction)
        % RV minor axis (in the anterior to posterior direction)
        % - Get the standard deviation of the blood pool in the LV to RV
        % direction, as a surrogate of the RV size:        
        RV_p = getPoints_fromBinary(RV,p0,[dx dx dx]);
        RVlength = GetStandardDeviationInDirection(RV_p,FirstAxis);
        RVmajorAxis = GetStandardDeviationInDirection(RV_p,SecondAxis);
        RVminorAxis = GetStandardDeviationInDirection(RV_p,ThirdAxis);
        
        % Now get the dimensions required for the synthesis of the
        % templates:
        dimensions.RVsizeminor = 1 + RVminorAxis/LVaxis;
        dimensions.RVsizeM     = 1 + RVmajorAxis/LVaxis;
        dimensions.RVapex      = 1 - RVlength/LVlength;

        
        % Flat to indiciate whether the RV apex is wide or slim:
        RVapexWidth = GetApicalBasalWidths(RV_p,R);
        LVapexWidth = GetApicalBasalWidths(LV_p,R);
        if RVapexWidth/LVapexWidth > 1
            dimensions.bBiVw = 1;
        else
            dimensions.bBiVw = 0;
        end
        % Quality control:
        if dimensions.RVapex < 0 || dimensions.RVapex > 1
            fprintf('ERROR! Analysis of the binary mask led to a location of the RV apex out of range (RVapex = %1.2f)!!\n',dimensions.RVapex);
            dimensions.RVapex = 0.5;
        end
    end    
end


end

function [aW,bW] = GetApicalBasalWidths(points,R)
% Function to get the width of the apical region of a cloud of points,
% oriented as the matrix R describes:

    [nPoints dim] = size(points);
    if nPoints<dim
        points = points';
        [nPoints dim] = size(points);
    end
    
    % 1. Define the apical points: those < mean-std of the main axis:
    FirstAxis = R(:,1);
    [S Coordinates] = GetStandardDeviationInDirection(points,FirstAxis);
    MeanC = mean(Coordinates,1);
    I1 = find(Coordinates < MeanC - S);
    I2 = find(Coordinates > MeanC + S);
    
    Points1 = points(I1,:);
    Points2 = points(I2,:);
    
    % 2. Find out the width fo the points in the two perpendicular axis to
    % the main axis:
    S2 = GetStandardDeviationInDirection(Points1,R(:,2));
    S3 = GetStandardDeviationInDirection(Points1,R(:,3));    
    W1 = sqrt(S2*S3);   
    S2 = GetStandardDeviationInDirection(Points2,R(:,2));
    S3 = GetStandardDeviationInDirection(Points2,R(:,3));    
    W2 = sqrt(S2*S3); 
    if W1<W2
        aW = W1;
        bW = W2;
    else
        aW = W2;
        bW = W1;
    end
end

function [S Coordinates] = GetStandardDeviationInDirection(points,direction)
    [nPoints dim] = size(points);
    if nPoints<dim
        points = points';
        [nPoints dim] = size(points);
    end
% It's only the std, no need to remove the mean!    
%     MeanPoint = mean(points,1);
%     Points = points - repmat(MeanPoint,nPoints,1);
    Coordinates = points * direction;
    S = std(Coordinates);
end

function [R,S] = Template_Initial_Alignment(seg_p,params)
% Due to the nature of how heart templates are synthetised, these are
% always aligned in the Z axis (up), and with the RV towards the Y axis:

    bOffsetXangleinBiValignment = params.bOffsetXangleinBiValignment;
    OffsetXangle = params.OffsetXangle;
    if bOffsetXangleinBiValignment
        % v1, aligned with LV main axis, is tilted towads the free wall:
        v12 = -sin(pi*OffsetXangle/180);
        v13 = cos(pi*OffsetXangle/180); 
    else
        v12 = 0; v13 = 1;        
    end
    v22 = v13;
    v23 = -v12;
    v1 = [0 v12 v13];
    v2 = [0 v22 v23];
    PCA= PCAtransformation();
    R = PCA.BuildRotMatrixFromTwoVectors(v1,v2);
    S = GetScaleFromPointsAndAxis(seg_p,R);
end

function v = normalise(v)
    v = v /  sqrt(sum(v.^2));
end

function [R,S] = Assemble_Horizontal_Initial_Alignment(seg_p,RVdir)
    % The main inertial axis is the z axis:
    v1 = [0 0 1];
    % The the other two axis are defined by RVdir:
    RVdir = RVdir / sqrt(dot(RVdir,RVdir));
    v2 = [RVdir 0];
    v3 = cross(v1,v2);
    R = [v1' v2' v3'];
    S = GetScaleFromPointsAndAxis(seg_p,R);
end
function S = GetScaleFromPointsAndAxis(seg_p,R)
    S = zeros(3,3);
    [nPoints b] = size(seg_p);
    for iCoor=1:3
        vect = R(:,iCoor);
        VectMatrix = repmat(vect',nPoints,1);
        coor = dot(VectMatrix',seg_p');
        S(iCoor,iCoor) = var(coor);
    end
end

function [R,S] = Horizontal_Inital_Alignment(seg_p)
    vx = [1 0 0];
    vy = [0 1 0];
    vz = [0 0 1];
    % The scale is obtained by the variance in each coordinate, since the
    % orientation is the coordinate axis:
    S = zeros(3,3);
    variance = zeros(1,3);
    for iCoor=1:3
        variance(iCoor) = var(seg_p(:,iCoor));
    end
    % It should be that in Z axis the variance is the biggest,
    % accordintly to the supposed main axis:
    if ~((variance(3)>variance(2)) && variance(3)>variance(1))
        fprintf('WARNING! The shape with low resolution in Z axis did not have the biggest variance in this axis!\n')
    end
    S(1,1) = variance(3);
    %Then, the second vector is the one with biggest variance:
    if(variance(2)>variance(1))
        S(2,2) = variance(2);
        S(3,3) = variance(1);
        R = [vz' vy' -vx'];
    else
        S(2,2) = variance(1);
        S(3,3) = variance(2);
        R = [vz' vx' vy'];
    end
end

function seg_im = find_myo(seg_im,seg_hd,seg_name,debug)

segs = unique(seg_im);
n = numel(segs);
disp('')
disp([seg_name ' has ' num2str(n-1) ' labeled segments']);
disp('')

% Distinguish myocardial segments
mass = zeros(length(segs)-1,1);av_d = mass;
tot = numel(find(seg_im>0));

for i = 1:length(segs)-1
    % based on mass
    mass(i) = numel(find(seg_im==(segs(i+1)))); % consider only non-zero segments
    im_tmp = seg_im==segs(i+1);
    
    % Based on mean dist.
    if mass(i)/tot < 0.1    % discard segments <10% myocardial mass
        av_d(i) = 0;
    else
        p_tmp = getPoints_fromBinary(im_tmp,seg_hd.origin,seg_hd.spacing);
        av_d(i) = mean_dist_from_COM(p_tmp);
    end
end

seg_im = seg_im == segs(find(av_d==max(av_d))+1);

if debug
    %imagesc(seg_im(:,:,round(size(seg_im,3)/2)));
    show_segment_surface(seg_im,seg_hd.origin,seg_hd.spacing,0.5,0.7);
    input('Press enter to close figure and continue in debug mode.');
    close;
end

end

function [R,S,bSwapLV_RV] = get_rotation_matrix(seg_type, LV, p0, dx, Seg_p, options)
bDebug = 0;
bRefRVdirection = 0;
bRVpoolLabel = 0;
bSwapLV_RV = 0;
% Default options
bZMainAxis = 0;
bBasalPoint = 0;

if nargin==6
    if isfield(options,'LVbasalPoint'), 
        BasalPoint = options.LVbasalPoint;
        bBasalPoint = 1;
    end
    if isfield(options,'RVbasalPoint'), 
        BasalPoint = options.RVbasalPoint;
        bBasalPoint = 1;
    end
    if isfield(options,'RVdirection'),
        RVdirection = options.RVdirection; 
        bRefRVdirection = 1;
    end
    if isfield(options,'RVbp'),
        % This is the binary mask of the RV blood pool
        RVbp = options.RVbp;
        bRVpoolLabel = 1;
    end
    if isfield(options,'bZMainAxis'), bZMainAxis = options.bZMainAxis; end
    if isfield(options,'RV'), RV = options.RV; end
    if isfield(options,'LVdir'), LVdir = options.LVdir; end        
end

% Check no contradictory options:
if exist('LVdir','var') && bZMainAxis
    fprintf('WARNING! contradictory parameters in definition of the main axis of inertia of the heart.\n');
    fprintf('         introduced axis of inertia prevaleces over assumption of axis = Z axis\n');
end
obj = PCAtransformation;
LV_p = getPoints_fromBinary(LV,p0,[dx dx dx]);
[foo, LV_V, LV_C] = obj.PrincipalComponentAnalysisOfShape(LV_p);

% There are two possibilities to overwrite the LV orientation:
% - Assumption that LV main axis is the Z axis (flag bZMainAxis)
% - LV main axis provided in LVdir
if(bZMainAxis)
    % Do not make a PCA, which might lead to a very tiltled main axis of
    % the LV if the shape is nearly spherical and has been truncated:
    vz_temp = [0 0 1]';
    vy_temp = [0 1 0]';
    vx_temp = [1 0 0]';
    LV_V = [vz_temp vy_temp -vx_temp];
    [LV_C] = obj.meanCoordinate(LV_p);
else    
    vz_temp = LV_V(:,1);             % define z-axis: PC1 of LVBP
end
if exist('LVdir','var')
    vz_temp = LVdir;    
    % Find any other two orthonormal vectors to this main axis (irrelevant
    % for the rest of the analysis, but at least a valid orientation
    LV_V = BuildOrientationMatrixFromMainVector(vz_temp);
end


% Ensure PC1 points toward base
t_pts = transform_coords_z_align(Seg_p,LV_V,LV_C);
if(bBasalPoint)
    centre = mean(t_pts,1);
    AproxVZ = BasalPoint - centre;
    if (dot(AproxVZ,vz_temp)<0)
        vz = - vz_temp;
    else
        vz = vz_temp;
    end
else
    [r] = transform_coords_cylindrical(t_pts);
    vz = check_PC1_direction(r, t_pts, vz_temp);
    % hack when dilated apex:
    % vz = vz_temp;
end

% define y-axis: use COMs of RV and LV
if strcmp(seg_type,'BiV') || bRVpoolLabel    
    if(bRVpoolLabel)
        fprintf('Computing orientation from the label given to the RV blood pool\n');
        RV_p = getPoints_fromBinary(RVbp,p0,[dx dx dx]);
    else
        RV_p = getPoints_fromBinary(RV,p0,[dx dx dx]);
    end
    [RV_C] = obj.meanCoordinate(RV_p);
    vy = (RV_C - LV_C)';
    if(bRefRVdirection)
        % Flip vy (can be wrong if the RV has more voxels) if opposed to
        % the reference direction introduced (from valve planes, more
        % reliable):
        RefDir = RVdirection;
        D = dot(vy,RefDir);
        if D<0
            vy = -vy;
            % RV is the LV in reality:
            bSwapLV_RV = 1;
        end
    end
    vy = cross(cross(vz,vy),vz);
    vy = vy/norm(vy);       % define y-axis
    fprintf('Orientation LV to RV centre (perpendicular to main axis) is: %1.2f, %1.2f, %1.2f (after normalization).\n',vy);
    vx = cross(vy,vz);      % define x-axis
    
elseif strcmp(seg_type,'LV')
    vy = LV_V(:,2);
    vx = LV_V(:,3);
    if vz ~= vz_temp
        vx = -vx;
    end
end

% Rotation matrix from template space to segmentation space
vx = reshape(vx,3,1);
vy = reshape(vy,3,1);
vz = reshape(vz,3,1);
%R = [vx vy vz];
% CORRECTION BY P.LAMATA: THE EIGENVECTORS SHOULD BE ALIGNED IN IMPORTANCE,
% and the main axis is the one in Z, then the one in Y, and then the one in
% X:
R = [vz vy vx];
R = obj.SetRightHandRule(R);

% Obtain scaling from entire segmentation - Matt edit: 27/02/12
% S = obj.PrincipalComponentAnalysisOfShape(Seg_p);
% Obtain scaling with current orientation, only once the orientation is 
% corrected by the search of the RV and the constraint of an horizontal 
% base! - Pablo edit: 01/03/12
% Correction! The eigenvectors are stored per row, but the rotation matrix
% has the axis per column!!! So the scale is not correct using this:
% S = obj.GetScaleGivenOrientation(Seg_p,R);
% It needs:
S = GetScaleFromPointsAndAxis(Seg_p,R);

% Further correction, 16/02/2012: correct the scales in terms of these new main axes!
% Remove of this correction, 01/03/2012 (P.Lamata):
% Since the shape might not be oriented to its PCA, these scale numbers are
% simply aligned with the R vectors, and there's no need to short them
% again!:
% Svalues = sort(unique(S),'descend');
% for iC=1:3
%     S(iC,iC) = Svalues(iC);
% end

if(bDebug)
    figure;
    plot_new_axes(LV,p0,dx,S,LV_C,R,Seg_p);
    title(sprintf('Scales: %1.1f, %1.1f, %1.1f',sqrt(S(1,1)),sqrt(S(2,2)),sqrt(S(3,3))));
end

end

function plot_new_axes(LV,p0,dx,S,C,V,points)

%s = (mean(LV_S(LV_S~=0)));

x_axis = [C ; (C + reshape(V(:,1),1,3)*sqrt(S(1,1)))];
y_axis = [C ; (C + reshape(V(:,2),1,3)*sqrt(S(2,2)))];
z_axis = [C ; (C + reshape(V(:,3),1,3)*sqrt(S(3,3)))];

show_segment_surface(LV,p0,[dx dx dx],0.9,0.5);
hold on, plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'linewidth',2,'color','r');
hold on, plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'linewidth',2,'color','g');
hold on, plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'linewidth',2,'color','b');
subsample = ceil(numel(points)/5000);
plot3(points(1:subsample:end,1),points(1:subsample:end,2),points(1:subsample:end,3),'.');
end

function t_pts = transform_coords_z_align(pts,V,C)
% Transform coordinates to be aligned in z with the long axis of the LVBP

% Use PC1 as z-axis.
R = [V(:,3) V(:,2) V(:,1)];
t_pts = pts;
for i = 1:length(pts)
    t_pts(i,:) = pts(i,:)-C;        % Centre the point coordinates about COM
    t_pts(i,:) = (R'*t_pts(i,:)')';  % Rotate points so z is aligned with LVBP LA
end
end

function [r,theta] = transform_coords_cylindrical(pts)
% Transform points to cylinderical coordinates

theta = zeros(length(pts),1); r = theta;

for i=1:length(pts)
    r(i) = norm(pts(i,1:2));
    theta(i)=atan2(pts(i,2), pts(i,1));
    % Theta not used here
end
end

function PC1 = check_PC1_direction(r,t_pts,PC1)
% Check PC1 points toward base
%  Criterion: Basal points have a greater mean radius than apical points
r_above_COM = r(t_pts(:,3)>0);
r_below_COM = r(t_pts(:,3)<0);
if mean(r_above_COM) < mean(r_below_COM)
    PC1 = -PC1;
end
end

function [dist] = mean_dist_from_COM(pts)
% Calculates the mean cuadratic distance (in order to save time)
[nPts foo] = size(pts);
COM  =  mean(pts,1);
DIS  = pts-ones(nPts,1)*COM;
DIS  = sum(DIS.^2,1);
dist = mean(DIS);
end