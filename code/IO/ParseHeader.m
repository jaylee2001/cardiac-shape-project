function [hPH,Mv2w] = ParseHeader(info,options)
% Function to translate the info header retrieved by the ReadData3D to the
% header style adopted from the read_image_vtk2 inhereted from D.Barber.
% 
% OUTPUT:
% - hPH: header, with the fields
%       origin
%       spacing
%       dim
%       TransformMatrix: 3x3
%       Mv2w: 4x4 matrix (voxel to world coordinate transformation)
%       File: name of the file
%       Extension: to indicate the format
% - Mv2w: matrix of transformation from voxel to world coordinates
%
% Source of ReadData3D:
% http://www.mathworks.com/matlabcentral/fileexchange/28139-read-medical-data-3d
%

% A reference about the Mv2w matrix: :http://www.mevislab.de/docs/2.1/MeVisLab/Resources/Documentation/Publish/SDK/GettingStarted/ch11s04.html 
% The default is the identity matrix:
Mv2w=zeros(3,4);
bRebuildMv2w = 1;
bRebuildTransfMat=0;
if nargin<2
    bForceRecomputeMatrices = 0;
    bNotReloadHeader = 0;
else
    if isfield(options,'bForceRecomputeMatrices'),  bForceRecomputeMatrices = options.bForceRecomputeMatrices; end
    if isfield(options,'bNotReloadHeader'),  bNotReloadHeader = options.bNotReloadHeader; end
end
for i=1:3
    Mv2w(i,i)=1;
end

%% Identify where data is coming from
loader = NaN;
if isfield(info,'axis'),        loader = 'scinrrd_load';        end
if isfield(info,'head'),        loader = 'read_image_vtk2';     end
if isfield(info,'hdr'),         loader = 'load_untouch_nii';    end
if isfield(info,'Filename'),    loader = 'ReadData3D';          end
if isfield(info,'ics_version'), loader = 'ICS';                 end
% Very raw formats:
FNs = fieldnames(info);
if length(FNs)==3
    if isfield(info,'origin')&&isfield(info,'spacing')&&isfield(info,'dim')
        loader = 'raw_format';
    end
end

%% Actions
switch loader
    case 'scinrrd_load'
        % Header coming from scinrrd_load
        hPH.spacing = [info(1).spacing info(2).spacing info(3).spacing];
        hPH.origin  = [info(1).min info(2).min info(3).min];
        hPH.dim     = [info(1).size info(2).size info(3).size];
        hPH.TransformMatrix = eye(3);
        hPH.CenterOfRotation= [0 0 0];    
    case 'raw_format'
        hPH = info;
        hPH.TransformMatrix = eye(3);
        hPH.CenterOfRotation= [0 0 0];
    case 'read_image_vtk2'
        % Header coming from read_image_vtk2: need to add transformation matrix
        fprintf(1,' ... parsing header of image %s (opened with Read_image_vtk2) ...',info.File);
        hPH = info;
        hPH.TransformMatrix = eye(3);
        hPH.CenterOfRotation= [0 0 0]; 
    case 'ICS'
        % This is an ICS header:
        hPH.spacing = info.parameter.scale(2:4);
        hPH.origin  = [0 0 0];
        hPH.dim     = info.layout.sizes(2:4);
        hPH.TransformMatrix = eye(3);
        hPH.CenterOfRotation= [0 0 0];
    case 'load_untouch_nii'
        % Header coming from load_untouch_nii
        fprintf(1,' ... parsing header of image %s (opened with load_untouch_nii) ...',info.fileprefix);
        hPH = info;
        if (~bNotReloadHeader)
            hPH.File = [info.fileprefix '.nii'];
            origin = [info.hdr.hist.qoffset_x info.hdr.hist.qoffset_y info.hdr.hist.qoffset_z];
            hPH.origin   = reshape(origin,1,3);
            hPH.spacing  = reshape(info.hdr.dime.pixdim(2:4),1,3);  
            hPH.dim      = reshape(info.hdr.dime.dim(2:4),1,3);
        
            % Now get the Mv2w matrix:
            NIIfile = [info.fileprefix '.nii'];
            hPH.TransformMatrix = eye(3);
            hPH.CenterOfRotation= [0 0 0];
            [Mv2w,TransformMatrix] = GetMv2wMatrixFromNifti(NIIfile);       
            if ~isnan(Mv2w(1)),   bRebuildMv2w = 0; hPH.Mv2w = Mv2w;    end
            if ~isnan(TransformMatrix(1)),   bRebuildTransfMat = 0; hPH.TransformMatrix = TransformMatrix;   end    
        end
        
    case 'ReadData3D'
        % Header coming from ReadData3D:
        fprintf(1,' ... parsing header of image %s (opened with ReadData3D) ... ',info.Filename);
        if ~isfield(info,'Format')
            % NIFTI header by ReadData3D does not return any Format field...
            if isfield(info,'Filename')
                [foo info.Format] = RemoveExtensionFromImageName(info.Filename);
            end
        end

        switch info.Format
            case {'gipl','GIPL','.gipl','.GIPL'}      
                hPH         = info;
                if isfield(info,'Origin')
                    hPH.origin  = reshape(info.Origin,1,3);
                else if isfield(info,'Origing')
                        % A hack to recover from an apparent bug in ReadData3D
                        hPH.origin  = reshape(info.Origing,1,3);
                    else 
                        fprintf('ERROR! No origin found in header\n')
                        hPH.origin = [0 0 0];
                    end
                end
                hPH.spacing = reshape(info.PixelDimensions,1,3);
                hPH.dim     = reshape(info.Dimensions,1,3);
                hPH.TransformMatrix = eye(3);
                fprintf('WARNING! Rotation matrix from GIPL file is lost! (need to code it in ParseHeader.m\n');
                % More info about this orientation in: http://www.cs.ucl.ac.uk/staff/X.Zhuang/zxhproj/classzxh_image_gipl_t.html#_details
                % image reader read gipl image, need to swap according to magic number Orientation issue: read: will consider header[106,110,114] and header[122,126,130] as x-axis and y-axis if header[138,142,146]==0 (method2), else header[106-154] is image to world matrix(method3) write: always use method2 unless det(imagetoworldmatrix)<0 then use method3 ADD: read analyze format (.hdr,.img),y-dim index should be mirror flip in order to match gipl
                %EXAMPLE:
    %              Filename: 'E:\HFfromManav\Case6/manual_lv.gipl'
    %             FileSize: 370556
    %           Dimensions: [115 115 14]
    %      PixelDimensions: [1.4015 1.4015 10.0000]
    %            ImageType: 16
    %              Patient: [1x80 char]
    %               Matrix: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    %          Orientation: 0
    %             VoxelMin: 0
    %             VoxelMax: 0
    %               Origin: [-57.8167 99.7146 -51.7006]
    %     RescaleIntercept: 0
    %         RescaleSlope: 0
    %        InterSliceGap: 0
    %             UserDef2: 0
    %                 Par2: 0
    %               Offset: 256
    %               Format: 'gipl'
            case {'.nii','.NII'}            
                hPH          = info;
                hPH.origin   = reshape([info.QoffsetX info.QoffsetY info.QoffsetZ],1,3);
                hPH.spacing  = reshape(info.PixelDimensions(1:3),1,3);   
                hPH.dim      = reshape(info.Dimensions(1:3),1,3);  
                hPH.TransformMatrix = eye(3);
                fprintf('WARNING! Rotation matrix from NIFTI file is lost! (need to code it in ParseHeader.m\n');
                % Example: 
    %                        Filesize: 230848
    %            Filename: 'D:\EuHeart\data\cases\JiaheCases\28\manual_lv2.nii'
    %           SizeofHdr: 348
    %            DataType: 512
    %              DbName: '                  '
    %             Extents: 0
    %        SessionError: 0
    %             Regular: 'r'
    %             DimInfo: ' '
    %          Dimensions: [98 98 12 1 1 1 1]
    %         headerbswap: 0
    %            IntentP1: 0
    %            IntentP2: 0
    %            IntentP3: 0
    %          IntentCode: 0
    %         datatypestr: 'UNKNOWN'
    %            bitvoxel: 0
    %         DataTypeStr: 'UINT16'
    %            BitVoxel: 16
    %              Bitpix: 16
    %          SliceStart: 0
    %     PixelDimensions: [7x1 double]
    %           VoxOffset: 352
    %            SclSlope: 1
    %            SclInter: 0
    %            SliceEnd: 0
    %           SliceCode: ' '
    %           XyztUnits: 2
    %       xyzt_unitsstr: 'UNKNOWN'
    %        XyztUnitsStr: 'MM'
    %              CalMax: 0
    %              CalMin: 0
    %      Slice_duration: 0
    %             Toffset: 0
    %               Glmax: 0
    %               Glmin: 0
    %             Descrip: '                                                                                '
    %             AuxFile: '                        '
    %           QformCode: 2
    %           SformCode: 1
    %            QuaternB: 0
    %            QuaternC: 0
    %            QuaternD: 1
    %            QoffsetX: -6.7174
    %            QoffsetY: 134.1687
    %            QoffsetZ: -120.0330
    %               SrowX: [4x1 double]
    %               SrowY: [4x1 double]
    %               SrowZ: [4x1 double]
    %          IntentName: '                '
    %               Magic: 'n+1 '

            case {'vtk','VTK','.vtk','.VTK'}
                hPH.head     = [info.Format ' DataFile Version ' info.Version];
                hPH.contents = info.Header;
                hPH.datatype = info.DatasetFormat;
                hPH.dataset  = info.DatasetType;
                hPH.origin   = info.Origin;
                hPH.spacing  = info.PixelDimensions;
                hPH.dim      = info.Dimensions;
                hPH.points   = prod(hPH.dim);
                hPH.scalars_name = info.DataName;
                hPH.format   = info.DataType;
                hPH.comp     = info.NumberOfComponents;
                % A rotation information is added (as in a MHA image format):
                hPH.TransformMatrix = eye(3);
                hPH.CenterOfRotation= [0 0 0];        

            case {'MHA','mha','.MHA','.mha'}
                % In this case, just copy all fields, and create duplicate ones
                % with the same key names:
                hPH          = info;
                hPH.origin   = info.Offset;
                hPH.spacing  = info.PixelDimensions;   
                hPH.dim      = info.Dimensions;  
                % Build the transformation matrix in a matrix form:
                if length(info.TransformMatrix)==9
                    % This is a vectorial indexing coming from a medical image that should
                    % be transformed into a matrix
                    temp = info.TransformMatrix;
                    RotationMatrix(1,1:3) = temp(1:3);
                    RotationMatrix(2,1:3) = temp(4:6);
                    RotationMatrix(3,1:3) = temp(7:9);
                    % TODO: CHECK IF THIS IS THE RIGHT WAY TO REBUILD THE MATRIX, should it
                    % be instead by columns? 
                    % 08/09/10: this seems correct, and a rotation matrix should be
                    % symetric!               
                    % Important bug solved: the rotation matrix in this Mv2w should have
                    % been with the permuted! (22/01/2013)
                    hPH.TransformMatrix = RotationMatrix';
                end
                % Example of a header of a MHA image:
                %                  Filename: 'D:\EuHeart\data\HumanFibres\domain.mha'
                %                    Format: 'MHA'
                %            CompressedData: 'false'
                %                ObjectType: 'image'
                %        NumberOfDimensions: 3
                %                BinaryData: 'true'
                %                 ByteOrder: 'false'
                %           TransformMatrix: [0.7718 8.5167e-010 0.6359 0.3349 0.8501 -0.4065 -0.5405 0.5267 0.6561]
                %                    Offset: [37.6990 -110.8630 -33.9493]
                %          CenterOfRotation: [0 0 0]
                %     AnatomicalOrientation: 'RAI'
                %           PixelDimensions: [1.2061 1.2061 1.2061]
                %                Dimensions: [76 86 112]
                %                  DataType: 'ushort'
                %                  DataFile: 'LOCAL'
                %                  BitDepth: 16
                %                HeaderSize: 402
            otherwise
                fprintf(1,'ERROR in ParseHeader!! Image type %s not implemented yet.\n',info.Format);
        end
        hPH.File = info.Filename;
    otherwise
        fprintf('WARNING: no image loader recognised\n');
        hPH = info;
end

ScaleMatrix = zeros(3,3);
InvScaleMat = zeros(3,3);
for i=1:3
    ScaleMatrix(i,i) = hPH.spacing(i);
    InvScaleMat(i,i) = 1/hPH.spacing(i);
end
if(bRebuildMv2w)|| bForceRecomputeMatrices
    % 24/06/2013: order inversed in this convolution:
    Mv2w(1:3,1:3) = hPH.TransformMatrix*ScaleMatrix;   
    Mv2w(1:3,4)   = hPH.origin;
    hPH.Mv2w=Mv2w;
end

if(bRebuildTransfMat)
    % need to calculate TransformaMatrix:
    hPH.TransformMatrix = ScaleMatrix \ hPH.Mv2w(1:3,1:3); % This is a faster way of the inverse, = inv(ScaleMatrix) * Mv2w
end
    
% Check correctness of the Transformation Matrix:
RR = hPH.TransformMatrix;
epsilon = 1e-12;
if det(inv(RR) - RR') > epsilon
    fprintf('ERROR(ParseHeader)! Invalid rotation matrix!\n');
    
end

% Check that the spacing is positive (recover from possible bugs):
if hPH.spacing(1)<0 || hPH.spacing(2)<0 || hPH.spacing(3)<0
    fprintf('WARNING! A negative value of spacing found: %1.1f,%1.1f,%1.1f (and corrected)\n',hPH.spacing);
    hPH.spacing = abs(hPH.spacing);
end

% Build the matrix of world to voxel coordinates:
Matrix = hPH.TransformMatrix*ScaleMatrix;
Mw2v(1:3,1:3) = inv(Matrix);
Mw2v(1:3,4)   = - inv(Matrix) * hPH.origin';
hPH.Mw2v = Mw2v;

% Addition of a copy of spacing, for resampling purposes (used in
% registration of sparse cases for example):
hPH.OriginalSpacing = hPH.spacing;

% Make sure the key parameters are in a row:
hPH.origin   = reshape(hPH.origin,1,3);
hPH.spacing  = reshape(hPH.spacing,1,3);  
hPH.dim      = reshape(hPH.dim,1,3);

% And put also the extension:
if isfield(hPH,'File')
    [foo hPH.Extension] = RemoveExtensionFromImageName(hPH.File);
else
    hPH.Extension = [];
end
fprintf(1,' Finished! The image has %i,%i,%i voxels.\n',hPH.dim);

