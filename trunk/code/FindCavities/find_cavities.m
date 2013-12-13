% Iterate erosion on mask to obtain LV and RV cavities
% Written by Matthew Sinclair
% Created: 01/11/11
% This script is called by cav_initial_alignment.m
%
% NOTE: Assumption that LV cavity volume is larger than RV cavity
% volume in all BiV segmentations (search 'TODO').
% ------------------------------------------------------------------------
% INPUTS:
%   (1) convex hull volume binary
%   (2) input segmentation binary (rebinned to spacing of (1))
%   (3) segmentation type ('LV', or 'BiV')
%   (4) origin
%   (5) spacing
%   (6) integer option for debugging (0 for none, >0 for debug)
% ------------------------------------------------------------------------
% OUTPUTS:
%
%   (1) c1 = LV bloodpool
%   (2) c2 = RV bloodpool (when dealing with BiV segmentations)
% ------------------------------------------------------------------------
%
% Revision control:
% 27/02/12: Created more robust criterion for distinguishing LV and RV
% cavities, using shape eccentricity (M.Sinclair)


function [c1 c2] = find_cavities(seg_ch_im_filled,seg_im_rebin,seg_type,p0,dx, debug)

if nargin<6
    debug = 0;
end

if strcmp(seg_type,'BiV')
    N_segs = 2;
    fprintf('\n Finding BiV segmentation cavities... \n');
elseif strcmp(seg_type,'LV')
    N_segs = 1;
    fprintf('\n Finding LV segmentation cavity... \n');
end

% Subtract seg from mask
seg_ch_im_filled_tmp = seg_ch_im_filled;
se = strel('arbitrary',1,1);        % define erosion neighbourhood [1]
cavs = logical(uint8(seg_ch_im_filled_tmp - seg_im_rebin));
[cavs2 N] = bwlabeln(cavs);  % label connected regions
cavs2 = uint8(cavs2);

% Number of iteration steps for erode
iter = 0;

MaxIter = 15;
while (N~=N_segs || iter == 0) && iter <= MaxIter
    

    % Create distance map of points from boundary (1s on boundary)
    seg_ch_im_filled_tmp = bwdist(logical(cavs2),'chessboard');
    seg_ch_im_filled_tmp = bwdist(seg_ch_im_filled_tmp,'chessboard');
    
    % Perform erode on boundary
    seg_ch_im_filled_tmp = imerode(seg_ch_im_filled_tmp,se);
    seg_ch_im_filled_tmp = logical(uint8(seg_ch_im_filled_tmp));
    iter = iter + 1;
    fprintf('Iteration #%i to find to big cavities...\n',iter);   
    if (iter == MaxIter)
        debug = 1;
    end
    % Subtract seg from eroded mask
    cavs = logical(uint8(seg_ch_im_filled_tmp - seg_im_rebin));
    
    % Label segments
    [cavs2 N] = bwlabeln(cavs);
    cavs2 = uint8(cavs2);
    
    % 'Clean' - remove small regions, <10% of non-zero binary.
    N_tmp1 = N;
    for i = 1:N
        n = numel(find(cavs2==i));
        if n < numel(find(cavs2>0))/10
            cavs2(cavs2==i) = 0;
            N_tmp = N_tmp1 - 1;
            N_tmp1 = N_tmp;
        end
    end    
    N = N_tmp1;
    
    % Visualise detected cavities -------------------------------------
    if debug
        k = unique(cavs2);
        subplot(121)
            show_segment_surface(seg_im_rebin,p0,[dx dx dx]);
            title('Rebinned binary mask')
        subplot(122)
        for i = 1:N
            cavs2_view = cavs2==k(i+1);
            show_segment_surface(cavs2_view,p0,[dx dx dx],0.2,0.6);
            if i ~=N
                hold on;
            end
        end
        title(['Iterations for erode + clean = ' num2str(iter)]);
        input('Continue and close figure [press enter]:','s'); close;
        figure; imagesc(squeeze(max(cavs2,[],3)));
        title(['MIP in Z. Erosion iteration number = ' num2str(iter)]);
        input('Continue and close figure [press enter]:','s'); close;
    end
    % -----------------------------------------------------------------

    if N == N_segs
        if strcmp(seg_type,'BiV')
            [c1 c2 N] = identify_cavities(cavs2,seg_type,N);
        elseif strcmp(seg_type,'LV')
            [c1] = identify_cavities(cavs2,seg_type);
        end
        
        if iter == 15
            error('Error: Max number of iterations (15) reached, cavities not found! Try entering debug mode.')
        end
    end
    
end 
  
    % Visualise detected cavities ----------------------------------------
    if debug
        fprintf(' Points in mask = %7.0f \n', numel(find(seg_ch_im_filled==1)));
        fprintf(' Points in mask after erosion = %7.0f \n', numel(find(seg_ch_im_filled_tmp==1)));
        show_segment_surface(seg_ch_im_filled,p0,[dx dx dx],0.9,0.2);
        hold on;
        if exist('c1','var')
            show_segment_surface(c1,p0,[dx dx dx],0.1,0.6);
            legend('Convex hull','LV Cavity');
            if exist('c2','var')
                show_segment_surface(c2,p0,[dx dx dx],0.5,0.6);
                legend('Convex hull','LV Cavity','RV Cavity');
            end
        end
        title(['Final cavity volume. Iterations for erode + clean = ' num2str(iter)]);
    end
    % --------------------------------------------------------------------
end


function [c1 c2 N] = identify_cavities(cavs2,seg_type,N)
% This function correctly identifies LV and RV cavity volumes.
% It tests several shape characteristics, including:
% (i) mass
% (ii) eccentricity - added 27/02/12 (new core criterion)
%
% Output:
% (1) c1 = LV
% (2) c2 = RV (when working with BiV cases)
% (3) N = number of disconnected volumes detected
debug = 1;
k = unique(cavs2);

% For BiV segmentations
if strcmp(seg_type,'BiV')
    c1 = cavs2==k(2);
    c2 = cavs2==k(3);
    
    % Calculate eccentricity - assumption that LV cavity has
    % eccentricity closer to 1 than RV cavity.
    ecc1 = calculate_eccentricity(c1);
    ecc2 = calculate_eccentricity(c2);
    
    if ecc1>ecc2
        c1 = cavs2==k(3);
        c2 = cavs2==k(2);
    end
    
    % Check relative volumes of cavities
    p = numel(find(c1>0))/numel(find(c2>0));
    if p > 1;
        p = 1/p;
    else
        disp('Warning! Volume of LV smaller than RV.')
    end
    if p<0.1;
        % Check that 'RV' is at least 10% volume of LV, if not then
        % assume that one of the volumes detected contains both the LV
        % and RV and the other volume is not needed. Carry on with
        % erosion iterations.
        N = 3;          % TODO: more robust check...
    else
        fprintf('\n LV and RV cavities found. \n');
    end
    
    % For LV segmentations
elseif strcmp(seg_type,'LV')
    c1 = cavs2==k(2);
    fprintf('\n LV cavity found. \n');
    % TODO: implement more robust check for LV geometry
end

end

function ecc = calculate_eccentricity(im)
% Calculate eccentricity of minor axes from PCA:
% closer to 1 = more ellipsoidal
% closer to 0 = circular

obj = PCAtransformation;
p = getPoints_fromBinary(im,[0 0 0],[1 1 1]);
[S,V,C] = obj.PrincipalComponentAnalysisOfShape(p);

ecc = sqrt(1-S(3,3)/S(2,2));

end