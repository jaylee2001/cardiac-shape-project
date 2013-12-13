function [V] = Build_Isosurface(binary,originORM_v2w,spacing,options)
% Create an isosurface from a binary image.
%
% Version control
% - 23/01/2013: add the MatlabOffset

CoordinatesSystem = 'PhysicalCoordinatesNeglectingOrientation';
MatlabOffset = [1 1 1];

% In order to account to the X-Y Row-Column permutation inherent in matlab:
bPermuteXY=1;
if nargin==1
    originORM_v2w = [0 0 0];
    spacing = [1 1 1];
    fprintf('Warning! isosurface built in voxel space (origin at [0 0 0] and spacing of [1 1 1]\n');
end
bIsovalueGiven=0;
if numel(originORM_v2w)>3
    CoordinatesSystem = 'WorldCoordinates';
end

switch CoordinatesSystem
    case 'WorldCoordinates'
        % The second argument is supposed to be M_v2w, the matrix to transform
        % voxel to world coordinates.
        M_v2w = originORM_v2w;
    case 'PhysicalCoordinatesNeglectingOrientation'
        % The second argument is supposed to be origin
        origin = originORM_v2w;
        if (~nargin==3)
            fprintf(1,'Error, not enouth arguments to build the isusurface!');
            V=NaN;
            return;
        end    
        if numel(origin)~=3
            fprintf('Error, origin has not 3 values (the x, y and z coordinates of this point)\n')
            return;
        end
end

if nargin==4
    if isfield(options,'ss')
        % Make the isosurface with a subsampled version of the image
        ss = options.ss;
        binary = binary(1:ss:end,1:ss:end,1:ss:end);
        switch CoordinatesSystem
            case 'WorldCoordinates'
                scaleMatrix = ss*eye(3);
                M_v2w(1:3,1:3) = scaleMatrix*M_v2w(1:3,1:3);
            case 'PhysicalCoordinatesNeglectingOrientation'
                spacing= spacing*ss;
        end
    end
    if isfield(options,'isovalue')
        level = options.isovalue;
        bIsovalueGiven =1;
    end
end

M = max(binary(:));
m = min(binary(:));
if ~bIsovalueGiven
    level= (M-m)/2;
end

fprintf(1,'     Building Isosurface... \n');
fprintf(1,'            ... with level=%1.1f of an image with max=%i and min=%i. Nvox0=%i, Nvox1=%i\n',level,M,m,numel(find(binary==0)),numel(find(binary==1)));
fprintf(1,'            ... and in coordinates: %s\n',CoordinatesSystem);
% Generate a vector with the X, Y an Z coordinates for the isosurface
SizeBin=size(binary);
if numel(SizeBin)<3
    fprintf(1,'ERROR!! in Build_Isosurface, image is not 3D! (size=%i,%i,%i)\n',size(binary));
    return;
end
if m==M
    fprintf(1,'ERROR!! in Build_Isosurface, image has no information (max=%i and min=%i)\n',M,m);
    return;
end

%ATTENTION! There is a different way isosurface understands the X-Y
%coordinates!!! It internally permutes them (as compared to a visualization
%of the 3D points generated by the binary image)

if(bPermuteXY)
    XYswapedbinary = permute(binary,[2 1 3]);
    binaryForIsosurface=XYswapedbinary; 
else
    binaryForIsosurface=binary;
end

switch CoordinatesSystem
    case 'WorldCoordinates'
%         xVox=0.5:1:(SizeBin(1)-0.5);
%         yVox=0.5:1:(SizeBin(2)-0.5);
%         zVox=0.5:1:(SizeBin(3)-0.5);
% Equivalent to the use of MatlabOffset:
        xVox = 0:1:(SizeBin(1)-1);
        yVox = 0:1:(SizeBin(2)-1);
        zVox = 0:1:(SizeBin(3)-1);
        if(bPermuteXY)
            V = isosurface(xVox,yVox,zVox,binaryForIsosurface,level);
        else
            V = isosurface(yVox,xVox,zVox,binaryForIsosurface,level);
        end
        % Now, transform vertices from voxel to world coordinates:
        [s1,s2] = size(V.vertices);
        nVert   = numel(V.vertices)/3;
        VoxCoords  = reshape(V.vertices,nVert,3);
        %offs = repmat(MatlabOffset,nVert,1);
        %VoxCoords = VoxCoords - offs;
        VoxCoords(:,4)= ones(nVert,1);
        points  = VoxCoords * M_v2w';
        V.vertices=reshape(points,s1,s2);
    case 'PhysicalCoordinatesNeglectingOrientation'      
        %canvas = zeros(size(binary));
        %canvas(:,1,1)=1;
        %Xpoints = getPoints_fromBinary(canvas,origin,spacing);
        %Xcoords = Xpoints(:,1);
        Xcoords = origin(1) : spacing(1) : origin(1)+(spacing(1)*(SizeBin(1)-1));
        Ycoords = origin(2) : spacing(2) : origin(2)+(spacing(2)*(SizeBin(2)-1));
        Zcoords = origin(3) : spacing(3) : origin(3)+(spacing(3)*(SizeBin(3)-1));
        V = isosurface(Xcoords,Ycoords,Zcoords,binaryForIsosurface,level);
end




