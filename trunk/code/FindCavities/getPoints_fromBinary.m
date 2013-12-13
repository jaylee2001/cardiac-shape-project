function [points] = getPoints_fromBinary(binary,originORM_v2w,spacing,RotationMatrix)
% Function to convert a binary image into a cloud of 3D points.
% Parameters, two options:
% - binary: 3D binary image of the object (any value > 0 is considered 1)
% OPTION 1:
% - origin: 3D coordinate of the initial point for sampling
% - spacing: resolution in x,y,z used for the sampling
% - RotationMatrix: indicates the "direction in which the
% voxels are counted in each of the three dimensions", reflecting the image
% orientation (RAS, etc...)
% OPTION 2: 
% - M_v2w: the matrix of voxel to world coordinates transformation
%
% By Pablo Lamata, Oxford, June 2009
% Version control
% - 23/01/2013: change to a matrix multiplication, because the 1/2 voxel
% seems to be a bug!
% - 10/05/2012: Add the half a voxel location!!

bGetVoxelCoordinatesAlso=0;
bChange2Matrix = 1;
MatlabOffset = [1 1 1];

ImplementationOption=1;
if nargin==1
    % Update on the 1st May 2012: the origin is 
    % -[0 0 0], in order to be consistent with Build_Isosurface
    % - Not [0.5 0.5 0.5], the centre of the first voxel, 
    % - Not [1 1 1], as it was before. The reason was the difference
    % between the location of the isosurface and the cloud of points ot
    % calculate the fitting accuracy observed in the JR dataset:
    origin = [0 0 0];% [1 1 1];
    spacing= [1 1 1];
    RotationMatrix=eye(3);
else
    if nargin==2
        ImplementationOption=2;
        M_v2w = originORM_v2w;
    else
        origin = originORM_v2w;
        if nargin<4
            RotationMatrix=eye(3);
        end
        if(bChange2Matrix)
            % make the same as if we had M_v2w:
            hd.spacing = spacing;
            hd.origin  = origin;
            hd.dim = size(binary);
            hd.TransformMatrix = RotationMatrix;
            hd.head = 'Fake header for Mv2w computation';
            hd.File = 'none';
            hd = ParseHeader(hd);
            M_v2w = hd.Mv2w;
        end
    end
end
fprintf(1,'        Starting conversion from binary to point coordinates... ');    
if (max(max(max(binary))))<=0 
    fprintf(1,'\n        WARNING! Image looks empty (no value greater than 0) ');    
end
Indexes = find (binary>0);
BinSize=size(binary);
n=0;
points=zeros(length(Indexes),3);
if(bGetVoxelCoordinatesAlso)
    Vpoints=zeros(length(Indexes),3);
end
for i=1:BinSize(1)
    for j=1:BinSize(2)
        for k=1:BinSize(3)
            if binary(i,j,k)>0
                n=n+1;
                switch ImplementationOption
                    case 1
                        if(bChange2Matrix)
                            bMatrixMultiplication = 1;                            
                        else
                            bMatrixMultiplication = 0;
                            p(1) = (i-1) * spacing(1);
                            p(2) = (j-1) * spacing(2);
                            p(3) = (k-1) * spacing(3);
                            points(n,1:3) = (origin) + [0.5 0.5 0.5].*spacing + p * RotationMatrix;
                        end
                    case 2
                        bMatrixMultiplication = 1;
                end
                if(bMatrixMultiplication)
                    VoxPoint = [i j k] - MatlabOffset;
                    VoxCoordinate = [VoxPoint 1];
                    bMatrixMultiplication = 1;
                    points(n,1:3)  = VoxCoordinate * M_v2w';
                    if(bGetVoxelCoordinatesAlso)
                        Vpoints(n,1:3) = VoxPoint;
                    end
                end
            end
        end
    end
end
if(bGetVoxelCoordinatesAlso)
    mean(Vpoints,1)
end
fprintf(1,'Conversion finished! %i points generated\n',n);