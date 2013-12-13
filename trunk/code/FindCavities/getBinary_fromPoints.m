function [BinaryImage,PointsOutOfRange] = getBinary_fromPoints(points,BinSize,origin,spacing)
% Function to generate a binary image from a cloud of 3D points. 
% Parameters:
% - points: 3D points of the object
% - BinSize: of the resulting image (number of voxels in x,y,z)
% - origin: 3D coordinate of the initial point for sampling
% - spacing: resolution in x,y,z used for the sampling
% By Pablo Lamata, Oxford, June 2009

fprintf(1,'        Starting conversion from point coordinates to binary... ');    
BinaryImage = false(BinSize);
PointsOutOfRange = 0;
TotalPoints = length(points);
for n=1:TotalPoints
    % TODO: check how to solve the 
    x = round((points(n,1)-origin(1))/spacing(1)) + 1;
    y = round((points(n,2)-origin(2))/spacing(2)) + 1;
    z = round((points(n,3)-origin(3))/spacing(3)) + 1;
    if ( x<1 || y<1 || z<1 || x > BinSize(1) || y > BinSize(2) || z > BinSize(3) )
        PointsOutOfRange = PointsOutOfRange + 1;
    else
        BinaryImage(x,y,z) = 1;
    end
end
fprintf(1,'Conversion finished.\n');
if PointsOutOfRange>0
    fprintf(1,'        WARNING!!!! IngetBinary_fromPoints a %2.2f%% of points were out of the volume!!!\n', 100*PointsOutOfRange/TotalPoints);
end

