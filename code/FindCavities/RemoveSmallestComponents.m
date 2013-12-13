function [im1,RemovalList] = RemoveSmallestComponents(image,MinNumPixels,connectivity,nCCpreserved)
    if nargin<2
        MinNumPixels = 10;
    end
    if nargin<3
        connectivity = 4;
    end
    if nargin<4
        nCCpreserved = 2;
    end
    im1 = image; 
    CC = bwconncomp(image,connectivity);
    if (CC.NumObjects>2)
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [foo,Index]=sort(numPixels,'descend');
        RemovalList = [];
        for i=nCCpreserved+1:CC.NumObjects
            iCC = Index(i);
            % Delete the voxels of the smallest groups:
            if(numPixels(iCC)<MinNumPixels)
                im1(CC.PixelIdxList{iCC}) = 0;
                RemovalList = [RemovalList CC.PixelIdxList{iCC}'];
            end
        end
    end