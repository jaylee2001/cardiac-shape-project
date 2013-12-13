function seg_ch_im_filled = BuildConvexHull(seg_p,seg_hd)
    % convex hull
    seg_ch_k = convhulln(seg_p);

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
end