function [V] = show_segment_surface(binary,originORM_v2w,spacing,c,alpha)
% Surface rendering of a cloud of points
% - seg: 3D binary mask
fprintf(1,'        Starting show_segment_surface... '); 
if nargin<2
    originORM_v2w = [0 0 0];
end
if nargin <3
    spacing = [1 1 1];
end
if nargin <4
    c = 0.1;
end
if nargin<5
    alpha = 0.3;
end

if isstruct(originORM_v2w)
    CoordinatesSystem = 'WorldCoordinates';
    M_v2w = originORM_v2w.Mv2w;
else
    if numel(originORM_v2w)>3
        CoordinatesSystem = 'WorldCoordinates';
        M_v2w = originORM_v2w;
    else
        CoordinatesSystem = 'PhysicalCoordinatesNeglectingOrientation';
    end
end

switch CoordinatesSystem
    case 'WorldCoordinates'
        % The second argument is supposed to be M_v2w, the matrix to transform
        % voxel to world coordinates.
        V = Build_Isosurface(binary,M_v2w);
    case 'PhysicalCoordinatesNeglectingOrientation'
        % The second argument is supposed to be origin
        origin = originORM_v2w;
        V = Build_Isosurface(binary,origin,spacing);
end

Show_Isosurface(V,c,alpha);
fprintf(1,'Finished!\n'); 




%% PREVIOUS VERSION FROM DAVID BARBER
% function mesh = show_segment_surface(seg,c,alpha,level)
% % Surface rendering of a binary mask
% % - seg: 3D binary mask
% 
% 
% if nargin < 4
%     level = 0.5;
% end
% V = isosurface(seg,level);
% mesh.nodes = V.vertices;
% mesh.elements = V.faces;
% mesh.scalars = ones(length(mesh.elements),1);
% if nargout == 1
%     return
% end
% w = mesh_to_iso(mesh);
% if nargin < 2
%     alpha = 1
%     c = 0.5;
% end
% if nargin < 3
%     alpha = 1;
% end
% w.FaceVertexCData = c;
% 
% %figure
% set(gcf,'Color','w');
% h = patch(w,'FaceColor','interp');view(3);camlight;lighting phong
% axis equal
% set(h,'EdgeColor','none')
% set(h,'FaceColor','flat')
% set(h,'FaceAlpha',alpha)
% set(h,'SpecularStrength',0.3)
% caxis([0 1])
% axis off
% 
% function V = mesh_to_iso(mesh)
% %
% V.vertices = mesh.nodes;
% V.faces = mesh.elements;