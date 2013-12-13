classdef PCAtransformation
    % TODO: check the right way to rotate the shapes!
   properties
       % NOTE: shape1 is the fixed one, shape2 is the transformed one.

        % - Centre1, Centre2: values to get the transformation and the points
        % around which make the rotations of the shapes.
        Centre1
        Centre2
        Tx
        % - V1, V2: Orientation matrices of each eigenanalysis. They will be needed
        % if the points want to be transformed after the PCA alignment (needed to
        % apply the right scaling). Rt = V1*V2';
        V1
        V2
        Rt 
        % Sc: Scale vector (3 scales corresponding to the directions of the
        % eigenvectors of the shapes).
        S1
        S2
        Sc
        
   end
   
   methods
       function obj = PCAtransformation()
           obj.Centre1=[0 0 0];
           obj.Centre2=[0 0 0];
           obj.V1=eye(3);
           obj.V2=eye(3);
           obj.Sc=[1 1 1];
           obj.Tx=[0 0 0];
           obj.Rt=eye(3);
       end
       function [obj,points1,points2] = Images_PCA_Registration(obj,binary1,binary2,header1,header2,options)
           % Affine registration between two binary images by a Principal Component
           % Analysis (PCA). The main assumption is that the shape of the object
           % represented in the binary images is similar. The static image is binary1
           % (the patient's heart)
           %
           % Input:
           % - binary1, binary2: 3D volumes of binary maps. Any value >0 is considered
           % to be 1. Each volume is simply a 3D matrix.
           % - header1, header2 [optional]: headers of the binary maps, with the
           % information about the origin and spacing of the images. 
           % - options(1): 
           %        1: (default) normal PCA
           %        2: alignment of centres of masses only
           %        3: alignment of centres of masses and orientation
           %        4: alignment of centres of masses and isotropic scale
           %        5: PCA alignment taking into accout the world
           %        coordinates of images (need the Mv2w in the header of
           %        both images)
           % [OPTIONAL]:
           % - options(2): bInvertInMajorAxis
           % - options(3): bInvertInSecondAxis
           % - options(4): bInvertInThirdAxis
           
           % Default options
                PCAopt  = [1 0 0 0];
                spacing1=[1 1 1];    origin1=[0 0 0];
                spacing2=[1 1 1];    origin2=[0 0 0];
                
           fprintf(1,'        Starting PCA_Registration... '); 
           if (nargin-1)<2
               error('Too few input arguments!');
           else
               if ((nargin-1)==2)
                   % Default options
               else
                   if ((nargin-1)==3)
                       spacing1=header1.spacing; origin1=header1.origin;
                   else
                       if ((nargin-1)==4)
                           spacing1=header1.spacing; origin1=header1.origin;
                           spacing2=header2.spacing; origin2=header2.origin;
                       else
                           if (nargin-1)==5
                               spacing1=header1.spacing; origin1=header1.origin;
                               spacing2=header2.spacing; origin2=header2.origin;
                               if numel(options)==4
                                   PCAopt = options;
                               else % fill with 0's
                                   nO = numel(options);
                                   PCAopt(1:nO) = options;
                               end
                           else
                               error('Too many input arguments!');
                           end
                       end
                   end
               end
           end
           if PCAopt(1)==5
               WichCoordinateSystem = 'WorldCoordinateSystem';
           else
               WichCoordinateSystem = 'PhysicalUnitsNeglectingOrientation';
           end
           % 1. Binary images are transformed to a cloud of points
           switch WichCoordinateSystem
               case 'WorldCoordinateSystem'
                    points1 = getPoints_fromBinary(binary1,header1.Mv2w);
                    points2 = getPoints_fromBinary(binary2,header2.Mv2w);
               case 'PhysicalUnitsNeglectingOrientation'
                    points1 = getPoints_fromBinary(binary1,origin1,spacing1);
                    points2 = getPoints_fromBinary(binary2,origin2,spacing2);
           end
           [obj] = obj.PCA_Registration_Mine(points1,points2,PCAopt);
       end
       
       function [obj] = PCA_Registration_Mine(obj,points1,points2,option)
           % 1. Mean difference between coordiantes: rigid translation
           %fprintf(1, '  Spacing of binary1 = (%d,%d,%d)\n  Spacing of binary2 = (%d,%d,%d)\n',spacing1(1),spacing1(2),spacing1(3),spacing2(1),spacing2(2),spacing2(3));
           %fprintf(1, '  Origin of binary1 = (%d,%d,%d)\n  Origin of binary2 = (%d,%d,%d)\n',origin1(1),origin1(2),origin1(3),origin2(1),origin2(2),origin2(3));
           %fprintf(1, '  Centre of binary1 = (%d,%d,%d)\n  Centre of binary2 = (%d,%d,%d)\n',Centre1(1),Centre1(2),Centre1(3),Centre2(1),Centre2(2),Centre2(3));

           [Sc1,Vec1,C1] = obj.PrincipalComponentAnalysisOfShape(points1);
           [Sc2,Vec2,C2] = obj.PrincipalComponentAnalysisOfShape(points2);
           obj.Centre1 = C1;
           obj.Centre2 = C2;
           obj.Tx = obj.Centre1 - obj.Centre2;

           obj.S1 = Sc1;
           obj.S2 = Sc2;
           obj.V1 = Vec1;
           obj.V2 = Vec2;
           obj.Rt = obj.V1*obj.V2';
           obj = obj.UpdateSc();
           obj    = obj.SetTrasnformationOption(option(1));
           if(option(2)), obj    = obj.InvertAxis(1);  end
           if(option(3)), obj    = obj.InvertAxis(2);  end
           if(option(4)), obj    = obj.InvertAxis(3);  end
           obj    = obj.SetTrasnformationOption(option(1));
       end
       function [obj] = Update(obj);
           obj = UpdateSc(obj);
           obj.Rt = obj.V1*obj.V2';
           obj.Tx = obj.Centre1 - obj.Centre2;           
       end
           
       function [obj] = UpdateSc(obj)
           obj.Sc = sqrt([obj.S1(1,1)/obj.S2(1,1) obj.S1(2,2)/obj.S2(2,2) obj.S1(3,3)/obj.S2(3,3)]);
       end
       function [obj] = InvertAxis(obj,whichAxis)
           % CONVENTION ADOPTED: flip the axis of first shape
           switch whichAxis
               case 1 %Invert axis 1, and 2 in order to maintain Right Hand Rule
                   obj.V1(:,1) = -obj.V1(:,1);
                   obj.V1(:,2) = -obj.V1(:,2);
               case 2 %Invert axis 2, and 3 in order to maintain Right Hand Rule
                   obj.V1(:,3) = -obj.V1(:,3);
                   obj.V1(:,2) = -obj.V1(:,2);
               case 3 %Invert only axis 3, no conservation of Right Hand Rule!
                   obj.V1(:,3) = -obj.V1(:,3);
               otherwise
                   fprintf('ERROR (in PCAtransformation)! Wrong selection of axis to invert\n');
                   pause
           end
                   
       end
       function [obj] = PCA_Registration_DowloadedCode(obj,points1,points2,option)
           % 1. Mean difference between coordiantes: rigid translation
           obj.Centre1 = obj.meanCoordinate(points1);
           obj.Centre2 = obj.meanCoordinate(points2);
           obj.Tx = obj.Centre1 - obj.Centre2;
           %fprintf(1, '  Spacing of binary1 = (%d,%d,%d)\n  Spacing of binary2 = (%d,%d,%d)\n',spacing1(1),spacing1(2),spacing1(3),spacing2(1),spacing2(2),spacing2(3));
           %fprintf(1, '  Origin of binary1 = (%d,%d,%d)\n  Origin of binary2 = (%d,%d,%d)\n',origin1(1),origin1(2),origin1(3),origin2(1),origin2(2),origin2(3));
           %fprintf(1, '  Centre of binary1 = (%d,%d,%d)\n  Centre of binary2 = (%d,%d,%d)\n',Centre1(1),Centre1(2),Centre1(3),Centre2(1),Centre2(2),Centre2(3));

           Points1WithoutMean = points1 - (ones(length(points1),1) * obj.Centre1);
           Points2WithoutMean = points2 - (ones(length(points2),1) * obj.Centre2);
           [U1,S1,V1] = pca(Points1WithoutMean,3);
           [U2,S2,V2] = pca(Points2WithoutMean,3);
                  %fprintf(1, '  Eigenvalues of binary1 = (%d,%d,%d)\n
                   %Eigenvalues of binary2 = (%d,%d,%d)\n',S1(1,1),S1(2,2),S1(3,3),S2(1,1),S2(2,2),S2(3,3));
           obj.S1 = S1;
           obj.S2 = S2;
           obj.V1 = V1;
           obj.V2 = V2;
           obj.Rt = obj.V1*obj.V2';
           obj.Sc = [sqrt(obj.S1(1,1)/obj.S2(1,1)) sqrt(obj.S1(2,2)/obj.S2(2,2)) sqrt(obj.S1(3,3)/obj.S2(3,3))];
           %obj.Sc = 1./obj.Sc;
           obj    = obj.SetTrasnformationOption(option);
       end
       
       function obj = SetTrasnformationOption(obj,option)
		   switch option
               case 1
                   % 3. Retrieve the transformation parameters for output
                   fprintf('        PCA alignment done complete. Scales = (%1.2f,%1.2f,%1.2f)\n',obj.Sc);
               case 2
                   obj.V1 = eye(3);
                   obj.V2 = eye(3);
                   obj.Rt = eye(3);   
                   obj.S1 = eye(3); 
                   obj.S2 = eye(3); 
                   obj.Sc = ones(1,3);
                   fprintf('        No PCA alignment done, only centres aligned. Translation = (%1.2f,%1.2f,%1.2f)\n',obj.Tx);
               case 3
                   obj.S1 = eye(3); 
                   obj.S2 = eye(3); 
                   obj.Sc = ones(1,3);
                   fprintf('        PCA alignment done complete. SCALES NEGLECTED\n');
               case 4
                   obj.V1 = eye(3);
                   obj.V2 = eye(3);
                   obj.Rt = eye(3); 
                   % Bugged version: scale should be averaged after sqrt!
                   %s1 = mean([obj.S1(1,1) obj.S1(2,2) obj.S1(3,3)]);
                   %s2 = mean([obj.S2(1,1) obj.S2(2,2) obj.S2(3,3)]);
                   s1 = sqrt(obj.S1);
                   s2 = sqrt(obj.S2);
                   meanS1 = mean([s1(1,1) s1(2,2) s1(3,3)]);
                   meanS2 = mean([s2(1,1) s2(2,2) s2(3,3)]);
                   obj.S1 = (meanS1^2)*eye(3); 
                   obj.S2 = (meanS2^2)*eye(3); 
                   obj = obj.UpdateSc();
                   fprintf('        NO PCA alignment done. Centres aligned and isotropic scale retrieved\n');
                   
           end
           %obj.V1 = obj.SetRightHandRule(obj.V1);
           %obj.V2 = obj.SetRightHandRule(obj.V2);
       end
       
       function [S,V,C] = PrincipalComponentAnalysisOfShape(obj,points)
           % Implementation following the "Computing PCA using the
           % covariance method" from
           % http://en.wikipedia.org/wiki/Principal_component_analysis
           bDebug=0;
           % 1. Calculate the empirical mean
           C = obj.meanCoordinate(points);
           % 2. Calculate the deviations from the mean
           N = length(points);
           if N<=3
               fprintf('Error in PrincipalComponentAnalysisOfShape, not enough points introduced!\n');
           end
           B = points - (ones(N,1) * C);
           % 3. Find the covariance matrix
           if (bDebug), fprintf('About to calculate covariance matrix...'); end
           COV = B'*B/N;
           if (bDebug), fprintf('Done! size = %i,%i\n',size(COV)); end
           if (bDebug), fprintf('About to calculate eigenvectors and eigenvalues...'); end
           [V,S] = eig(COV);
           % 4. Reorder vectors:
           ss = [S(1,1) S(2,2) S(3,3)];
           [a,I] = sort(ss,'descend');                
           Stemp = zeros(3,3);   Vtemp = zeros(3,3);
           for i=1:3
               Stemp(i,i) = S(I(i),I(i));
               Vtemp(:,i) = V(:,I(i));
           end
           S = Stemp;
           V = Vtemp;
           [V] = obj.SetRightHandRule(V);
           if (bDebug), 
               fprintf('Done!\n'); 
               V
           end
           
       end
       
       function S = GetScaleGivenOrientation(obj,points,V)
           bDebug = 0;
           % In some specific cases, the "scaling" of the shape to a given
           % set of eigenvectors is required:
           C = obj.meanCoordinate(points);
           N = length(points);
           if N<=3
               fprintf('Error in PrincipalComponentAnalysisOfShape, not enough points introduced!\n');
           end
           B = points - (ones(N,1) * C);   
           % Now calculate S as the variance of B in each of the directions
           % indicated by V:
           S = zeros(3,3);
           if(bDebug), figure; end
           for iC = 1:3
               % Get the eigenvector, each row of V:
               eigvec = V(iC,:);
               % Make the dot product of points in this direction:
               data = B * eigvec';
               % Calculate the variance in this direction:
               S(iC,iC) = var(data);
               if(bDebug), subplot(3,1,iC), hist(data); title(sprintf('Histogram data*eigenvector(%i)',iC)); end
           end
       end
           
       function [obj] = SetRotationByAngles(obj,RotX,RotY,RotZ,ConcatenationOrder)
           obj.V1 = CreateRotationMatrix(RotX,RotY,RotZ,ConcatenationOrder);
           obj.V2 = eye(3);
           obj.Rt = obj.V1*obj.V2';
       end
       function [obj] = SetV1V2(obj,V1,V2)
           if(obj.checkRotationMatrix(V1))
                obj.V1 = V1;
           else
               fprintf('ERROR! Not valid rotation matrix in V1\n');
           end
           if(obj.checkRotationMatrix(V2))
                obj.V2 = V2;
           else
               fprintf('ERROR! Not valid rotation matrix in V2\n');
           end
           obj.Rt = obj.V1*obj.V2';
       end
       function [obj] = SetV1byPoints(obj,Points)
           % Function not tested yet (03/12/2012)
           V1 = BuildOrientationMatrixByPoints(Points);
           obj.V1 = V1;
       end           
           
       function [obj] = SetRotationMatrixInV1(obj,RotMatrix)
           if(obj.checkRotationMatrix(RotMatrix))
                obj.V1 = RotMatrix;
                obj.V2 = eye(3);
                obj.Rt = obj.V1*obj.V2';
           else
               fprintf('Not able to set rotation matrix, it is not valid one!\n');
           end
       end
       function [obj] = SetRotationMatrixInV2(obj,RotMatrix)
           if(obj.checkRotationMatrix(RotMatrix))
               obj.V2 = RotMatrix;
               obj.V1 = eye(3);
               obj.Rt = obj.V1*obj.V2';
           else
               fprintf('Not able to set rotation matrix, it is not valid one!\n');
           end
       end
       function [obj] = SetS1(obj,S)
           obj.S1 = S;
           obj = obj.UpdateSc();
       end
       function [obj] = SetS2(obj,S)
           obj.S2 = S;
           obj = obj.UpdateSc();
       end
       function [obj] = SetScaleAnisotropyInS1(obj,ScaleMatrix)
           bMainVectorvsTheOtherTwo=0;
           anisotropy = 3*ScaleMatrix/sum(ScaleMatrix(:));
           %anisotropy = ScaleMatrix/norm(ScaleMatrix(:));
           if(bMainVectorvsTheOtherTwo)
               anisotropy(2:3) = ones(1,2)*mean(anisotropy(2:3));
           end
           obj.S1 = obj.S1.*anisotropy;
           obj = obj.UpdateSc();
       end
       function [obj] = SetScaleAnisotropyInS2(obj,ScaleMatrix)
           bMainVectorvsTheOtherTwo=0;
           anisotropy = 3*ScaleMatrix/sum(ScaleMatrix(:));
           %anisotropy = ScaleMatrix/norm(ScaleMatrix(:));
           if(bMainVectorvsTheOtherTwo)
               anisotropy(2:3) = ones(1,2)*mean(anisotropy(2:3));
           end
           obj.S2 = obj.S2.*anisotropy;
           obj = obj.UpdateSc();
       end
       function [obj] = RandomTinyModification(obj)
           % Make a random change of the transformation parameters:
           % 1. A random rotation:
           angles = rand(1,3);
           R = CreateRotationMatrix(angles(1),angles(2),angles(3));
           obj.V1 = obj.V1 * R;
           
           % 2. A random tiny scaling:
           scaling = [1 1 1] + rand(1,3)/1000;
           for iS = 1:3
               obj.S1(iS,iS) = obj.S1(iS,iS) * scaling(iS);
           end
           
           % Update:
           obj = obj.Update();
       end
       function [TransformedPoints] = transform(obj,points,SenseOption,DOFsOption)
           % Function to scale points in the directions of the PCA analysis
           % - points: 3D points. It should be (3 x npoints):
           % - option: 
           %     1: (default) Shape1 is fixed, and Shape2 the transformed one, as how the PCA_Registration is done
           %     2: the opposite.
           % - DOFsOption:
           %     1: (default) a complete transformation
           %     2: Translation and rotation only (no scale)
           %     3: Translation only (no rotation, no scale)
           %     4: Scale and rotation (no translation) - for derivatives
           %     in a Hermite mesh
           
           if nargin<3
               SenseOption=1;
           end
           if nargin<4
               DOFsOption=1;
           end
           switch SenseOption
               case 1
                   Cent2 = reshape(obj.Centre2,1,3);
                   Cent1 = reshape(obj.Centre1,1,3);
                   Scale = reshape(obj.Sc,1,3);
                   Vshp2 = obj.V2;
                   Vshp1 = obj.V1;
               case 2
                   Cent1 = reshape(obj.Centre2,1,3);
                   Cent2 = reshape(obj.Centre1,1,3);
                   Scale = [1 1 1]./(reshape(obj.Sc,1,3));
                   Vshp1 = obj.V2;
                   Vshp2 = obj.V1;
               otherwise
                   fprintf('Wrong selection of SenseOption in PCAtransformation.transform\n');
           end
           
           switch DOFsOption
               case 1
                   % Do nothing, defalut option
               case 2
                   % Turn scales to 1 (no scale wanted)
                   Scale = [1 1 1];
               case 3
                   % Turn scales to 1 (no scale wanted) and rotations to
                   % eye (no rotation wanted)
                   Scale = [1 1 1];
                   Vshp1 = eye(3);
                   Vshp2 = eye(3);
               case 4
                   % Turn Tx to 0
                   Cent2 = [0 0 0];
                   Cent1 = [0 0 0];
               otherwise
                   fprintf('Wrong selection of DOFsOption in PCAtransformation.transform\n');
           end
           nPoints=numel(points)/3;
           points2=reshape(points,nPoints,3);
           
           Points2WithoutMean = points2 - (ones(nPoints,1) * Cent2);
           %Changed by P.Lamata, 23/12/2011: remove the sqrt, Sc already has it!
           %MatrixOfScales     = ones(nPoints,1) * sqrt(Scale);
           MatrixOfScales     = ones(nPoints,1) * (Scale);
           Points2TranslatedRotatedScaled = MatrixOfScales .* (Vshp2'*Points2WithoutMean')';
           Points2TranslatedRotatedScaled = Vshp1*Points2TranslatedRotatedScaled';
           TransformedPoints  = Points2TranslatedRotatedScaled' + (ones(nPoints,1) * Cent1);
%            else
%                Points2WithoutMean = points2 - Cent2';
%                Points2TranslatedRotatedScaled = Scale'.*(Vshp2'*Points2WithoutMean);
%                Points2TranslatedRotatedScaled = Vshp1*Points2TranslatedRotatedScaled;
%                TransformedPoints = Points2TranslatedRotatedScaled + Cent1';
%            end
%            
% if nPoints>1
%     Points2WithoutMean = points2 - (ones(nPoints,1) * Centre2);
%     MatrixOfScales = ones(nPoints,1) * Sc;
%     Points2TranslatedRotatedScaled = MatrixOfScales .* (V2'*Points2WithoutMean')';
%     Points2TranslatedRotatedScaled = V1*Points2TranslatedRotatedScaled';
%     TransformedPoints = Points2TranslatedRotatedScaled' + (ones(nPoints,1) * Centre1);
% else
%     Points2WithoutMean = points2 - Centre2';
%     Points2TranslatedRotatedScaled = Sc'.*(V2'*Points2WithoutMean);
%     Points2TranslatedRotatedScaled = V1*Points2TranslatedRotatedScaled;
%     TransformedPoints = Points2TranslatedRotatedScaled + Centre1';
% end

           % NOTE: the need to scale at the "intermediate rotation" prevents this kind
           % of implementation:
           % Points2TranslatedRotatedScaled = MatrixOfScales .* (Rt*Points2WithoutMean')';
           % for i=1:length(points)
           %    TransformedPoints(i,:) =  Rt*Points2WithoutMean(i,:)'
       end
       
       function [map] = CreatePCAdeformationMap(obj,hImage)
           % Function to create a deformation map from the parameters of a PCA
           % alignment, in a set of coordinates defined in the header of an image.
           % INPUT:
           % - hImage:
           %     hImage.origin
           %     hImage.spacing
           %     hImage.dim 
           %--------------------------------
           % Parameters:
           % Subsampling of the image resolution:
           ss = min(floor(hImage.dim/50))+10;
           
           %--------------------------------
           nCoord = floor(hImage.dim/ss);
           % Make a +2 to avoid points out of the region of the DF
           ni=nCoord(1)+2;
           nj=nCoord(2)+2;
           nk=nCoord(3)+2;
           o = hImage.origin;
           s = hImage.spacing;
           map.x=zeros(ni,nj,nk);
           map.y=zeros(ni,nj,nk);
           map.z=zeros(ni,nj,nk);
           map.u=zeros(ni,nj,nk);
           map.v=zeros(ni,nj,nk);
           map.w=zeros(ni,nj,nk);
           fprintf('        Creating def field of size(%i,%i,%i) - subsampling factor %i - from PCA param: Sc(%1.1f,%1.1f,%1.1f), Tx(%1.1f,%1.1f,%1.1f)\n',ni,nj,nk,ss,obj.Sc,obj.Tx);
           
           % TODO: speed up making a matrix form of the loop
           %P1(1:ni,1:nj,1:nk) =  o + ss*[(1:ni)-1,(1:nj)-1,(1:nk)-1].*s;
           for i=1:ni
               for j=1:nj
                   for k=1:nk
                       P1= o + ss*[i-1,j-1,k-1].*s;
                       P2= obj.transform(P1);
                       map.x(i,j,k)=P1(1);
                       map.y(i,j,k)=P1(2);
                       map.z(i,j,k)=P1(3);
                       map.u(i,j,k)=P2(1)-P1(1);
                       map.v(i,j,k)=P2(2)-P1(2);
                       map.w(i,j,k)=P2(3)-P1(3);
                   end
               end
           end
           I = isnan(map.z);
           II = find(I);
           if numel(II)>0
               fprintf(1,'ERROR in PCAtransformation.CreatePCAdeformationMap!! Field has NaN!!\n');
           end
       end
       
       function [Sc1,Vec1,C1] = ShowMainAxisOfShape(obj,data,imageHeader,WhichCoords)
           % Function to show the main axis of a shape. Data is given by
           % either a list of points or by a binary image (and then the
           % third argument is needed)
           [a b c] = size(data);
           if nargin<4
               WhichCoords = 'PhysicalCoordinates_IgnoringOrientationOfImages';
           end
           if nargin<3
               DataGiven = 1; % cloud of points
           else
               DataGiven = 2; % an image
               % check that data is an image:
               if a==1 | b==1 | c==1
                   fprintf('Error in arguments of ShowMainAxisOfShape!\n');
               end
           end
           switch DataGiven
               case 2
                   % Need to convert image into cloud of points                   
                   switch WhichCoords
                       case 'PhysicalCoordinates_IgnoringOrientationOfImages'                           
                           origin1 = imageHeader.origin;
                           spacing1= imageHeader.spacing;
                           points1 = getPoints_fromBinary(data,origin1,spacing1);
                       case 'WorldCoordinates'
                           points1 = getPoints_fromBinary(data,imageHeader.Mv2w);
                   end
               case 1
                   points1 = data;
           end
           [Sc1,Vec1,C1] = obj.PrincipalComponentAnalysisOfShape(points1);
           figure
           hold on
           obj.draw_PrincipalAxes(C1,Vec1,Sc1,'g');
           switch DataGiven
               case 1
                   plot3(points1(:,1),points1(:,2),points1(:,3),'r.')
               case 2
                   % Create approximated isosurface (with subsampling)
                   options.ss = 3; % Subsampling
                   options.isovalue = 0.5; % Subsampling
                   switch WhichCoords
                       case 'WorldCoordinates'
                            NotNeeded=1;
                            V=Build_Isosurface(data,imageHeader.Mv2w,NotNeeded,options);
                       case 'PhysicalCoordinates_IgnoringOrientationOfImages'
                            V=Build_Isosurface(data,origin1,spacing1,options);
                   end
                   p=patch(V);
                   set(p,'FaceColor','red');%,'EdgeColor','none');
           end
           axis equal
       end
          
       function [] = ShowRegistrationResult(obj,points1,points2,senseOption)
           if nargin<4
               senseOption=1;
           end
            figure('Color',[1 1 1]) 
            subplot(221); hold on
            plot3(points1(:,1),points1(:,2),points1(:,3),'r.')
            obj.draw_PrincipalAxes(obj.Centre1,obj.V1,obj.S1,'b');
            plot3(points2(:,1),points2(:,2),points2(:,3),'g.')
            obj.draw_PrincipalAxes(obj.Centre2,obj.V2,obj.S2,'k');
            view(40,30); axis equal;
            title('Initial configuration');

            Points2Translated = obj.transform(points2,senseOption,3);
            subplot(222); hold on
            plot3(points1(:,1),points1(:,2),points1(:,3),'r.')
            obj.draw_PrincipalAxes(obj.Centre1,obj.V1,obj.S1,'b');
            plot3(Points2Translated(:,1),Points2Translated(:,2),Points2Translated(:,3),'g.')
            obj.draw_PrincipalAxes(obj.Centre2+obj.Tx,obj.V2,obj.S2,'k');
            view(40,30); axis equal;
            title('Rigid translation');

            Points2TranslatedRotated = obj.transform(points2,senseOption,2);
            Pts = Points2TranslatedRotated - (ones(length(points2),1) * obj.meanCoordinate(Points2TranslatedRotated));
            %[U1,S2r,V2r] = pca(Pts,3);
            [U1,S2r,V2r] = obj.PrincipalComponentAnalysisOfShape(Pts);
            subplot(223); hold on
            plot3(points1(:,1),points1(:,2),points1(:,3),'r.')
            obj.draw_PrincipalAxes(obj.Centre1,obj.V1,obj.S1,'b');
            plot3(Points2TranslatedRotated(:,1),Points2TranslatedRotated(:,2),Points2TranslatedRotated(:,3),'g.')
            obj.draw_PrincipalAxes(obj.Centre2+obj.Tx,V2r,S2r,'k');
            view(40,30); axis equal;
            title('Rigid translation and rotation');

            TransformedPoints2 = obj.transform(points2,senseOption,1);
            subplot(224); hold on
            plot3(points1(:,1),points1(:,2),points1(:,3),'r.')
            obj.draw_PrincipalAxes(obj.Centre1,obj.V1,obj.S1,'b');
            plot3(TransformedPoints2(:,1),TransformedPoints2(:,2),TransformedPoints2(:,3),'g.')
            obj.draw_PrincipalAxes(obj.Centre2+obj.Tx,V2r,S2r,'k');
            view(40,30); axis equal;
            title('Rigid translation, rotation and scaling');
       end
       function [InitialAlignmentTransformation] = TestFunction(obj)
           dir = 'D:\EuHeart\data\MedIA workbench\';
           %image1 = [dir 'plane1.vtk'];
           %image2 = [dir 'plane3.vtk'];
           image1 = [dir 'shape1.vtk'];
           image2 = [dir 'shape3.vtk'];
           image1 = [dir 'ss1.vtk'];
           image2 = [dir 'ss3.vtk'];
           image2 = 'D:\EuHeart\data\Jiahe\LVframe08cropped.vtk';
           image1 = 'D:\EuHeart\data\MedIA workbench\ellipsoids\LV1_binaryRes2.535_2.11_3.1448.vtk';
           image1 = 'D:\EuHeart\data\MICCAI challenge\Dataset_1-June28\1-Heart_Anatomy\2-Processed_Data\DT-MR_Data\Myocardium_Segmentation\img.vtk';
           image1 = 'D:\EuHeart\data\cases\case14\meshing\case14v1.vtk';
           image2 = 'D:\EuHeart\data\Heart models\Pablo\heartMedIAb_binaryRes2.1875_2.1875_3.6.vtk';
           image1 = 'N:\data\KCLcases\case15\meshing\case15Ventricles.vtk';
           image2 = 'D:\EuHeart\data\Heart models\Pablo\heartMedIAb_binaryRes2.1875_2.1875_3.6.vtk';
           
           [I1,hI1] = read_image_vtk2(image1); I1 = logical(I1);
           [I2,hI2] = read_image_vtk2(image2); I2 = logical(I2);
           InitialAlignmentTransformation = PCAtransformation();
           [InitialAlignmentTransformation,points1,points2] = InitialAlignmentTransformation.Images_PCA_Registration(I1,I2,hI1,hI2);
           InitialAlignmentTransformation.ShowRegistrationResult(points1,points2,1);
       end
       
       function V1 = BuildOrientationMatrixFromPoints(obj,Points)
           % Four points expected: first two points define the first axis
           % of intertia, the second two the second axis.
           if numel(Points)~=12
               fprintf('ERROR! the number of points introduced is not as expected!\n')
               fprintf('       numel(Points) = %i (and it should be 12)\n',numel(Points));
               fprintf('       in function "BuildOrientationMatrixByPoints" from PCAtransformation\n')
           else
               Pts = reshape(Points,4,3);
               v1 = Pts(2,:) - Pts(1,:);
               v2 = Pts(4,:) - Pts(3,:);
               V1 = obj.BuildRotMatrixFromTwoVectors(v1,v2);
           end
       end
       function R = BuildRotMatrixFromTwoVectors(obj,v1,v2)       
           v1 = obj.normalise(v1);
           v2 = obj.normalise(v2);
           v3 = cross(v1,v2);
           % For some reason, the module of these three vectors is not completely
           % 1:
           v3 = obj.normalise(v3);
           % And the initial v1 and v2 might not be completely
           % perpendicular:
           v2 = cross(v3,v1);
           v2 = obj.normalise(v2);
           R = [v1' v2' v3'];
           bCorrect = obj.checkRotationMatrix(R);
           if (~bCorrect)
               a=1;
           end
       end
   end
   
   methods (Static=true)
       function [V] = SetRightHandRule(V)
           % TODO: CHECK THAT THIS IS THE RIGHT IMPLEMENTATION. There must
           % be a more consistent way to check the correct orientation of
           % the eigenvectors for the rotation matrices!!! The objective
           % was to remove the "mirroring effect" that I was getting in the
           % second eigenvector (what happened in the heart and in the
           % plane shapes which I tried). But I am not sure at all of this.
           
           % Function to check that the principal vectors are aligned
           % following a right hand rule:
           v1 = V(:,1);
           v2 = V(:,2);
           v3 = V(:,3);
           
           v = cross(v1,v2);
           a = dot(v,v3);
           if a<0
               % The right hand rule is not respectd, need to invert v3
               % 05/06/2011(PL): inversion of V1 omitted!
               %V(:,1) = -v1;
               V(:,3) = -v3;
               fprintf('Right hand rule corrected!\n');
           end
       end
       function []  = draw_PrincipalAxes(centre,V,S,color)
            %coloraxes=color;
            % plot3(points(:,1),points(:,2),points(:,3),color)
            if nargin<3
                S = eye(3);
            end
            PrincipalAxe1(1,:) = centre; PrincipalAxe1(2,:) = centre + 3*V(:,1)' .* sqrt(S(1,1));
            PrincipalAxe2(1,:) = centre; PrincipalAxe2(2,:) = centre + 3*V(:,2)' .* sqrt(S(2,2));
            PrincipalAxe3(1,:) = centre; PrincipalAxe3(2,:) = centre + 3*V(:,3)' .* sqrt(S(3,3));

            plot3(PrincipalAxe1(:,1),PrincipalAxe1(:,2),PrincipalAxe1(:,3),'r','LineWidth',2);
            plot3(PrincipalAxe2(:,1),PrincipalAxe2(:,2),PrincipalAxe2(:,3),'g','LineWidth',2);
            plot3(PrincipalAxe3(:,1),PrincipalAxe3(:,2),PrincipalAxe3(:,3),'b','LineWidth',2);
       end
       function [centre] = meanCoordinate(points)
            centre(1)=mean(points(:,1));
            centre(2)=mean(points(:,2));
            centre(3)=mean(points(:,3));
       end
       function [v] = normalise(v)
           v = v/sqrt(sum(v.^2));
       end
       function [bCorrect] = checkRotationMatrix(R)
           epsilon = 10^(-5);
           if abs(det(R)-1)>epsilon
               fprintf('Wrong rotation matrix! Determinant = %1.3f\n',det(R));
               bCorrect = 0;
           else
               bCorrect = 1;
           end             
       end
       
   end
end
