function [img_siz] = io_GetNIIimageSize(hdr)

%  For each frame, precision of value will be read 
      %  in img_siz times, where img_siz is only the 
      %  dimension size of an image, not the byte storage
      %  size of an image.
      %
      img_siz = prod(hdr.dime.dim(2:8));

      %  For complex float32 or complex float64, voxel values
      %  include [real, imag]
      %
      if hdr.dime.datatype == 32 | hdr.dime.datatype == 1792
         img_siz = img_siz * 2;
      end
	 
      %MPH: For RGB24, voxel values include 3 separate color planes
      %
      if hdr.dime.datatype == 128 | hdr.dime.datatype == 511
	 img_siz = img_siz * 3;
      end