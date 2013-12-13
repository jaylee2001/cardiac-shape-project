function [] = Show_Isosurface(V,colour,alpha)
if nargin<2
    colour=1;
end
if nargin<3
    alpha = 1;
end
RGBcolor = ind2rgb(round(256*colour),jet(256));
h2 = patch(V,'FaceAlpha',alpha,'FaceColor',RGBcolor);
set(h2,'EdgeColor','none')
view(3);camlight;lighting phong
xlabel('X'),ylabel('Y'),zlabel('Z')
axis equal

%w=V;
%w.FaceVertexCData = colour;

%figure
%set(gcf,'Color','w');
% h = patch(w,'FaceColor','interp');view(3);camlight;lighting phong
% axis equal
% set(h,'EdgeColor','none')
% set(h,'FaceColor','flat')
% set(h,'FaceAlpha',alpha)
% set(h,'SpecularStrength',0.3)
% caxis([0 1])
% xlabel('X'),ylabel('Y'),zlabel('Z')