function [plane] = FitFlatPlane2Points(Points)
    [coeff,score,roots] = princomp(Points);
    basis = coeff(:,1:2);
    normal = coeff(:,3);
    plane(1:3)= normal;
    centre = mean(X,1);
    plane(4) = -plane(1:3)*centre';