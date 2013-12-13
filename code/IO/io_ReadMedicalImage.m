function [PH,hPH] = io_ReadMedicalImage(file,options)
% Function that encapsulates the reading of several medical imaging
% formats, and returns a simplified and unified one adapted to the use of
% the code by Pablo Lamata

bViewIsosurface = 0;

[nameHeart formatHeartFile] = RemoveExtensionFromImageName(file);
switch formatHeartFile
    case {'.bm'}
        ImReader = 'PhilipsImage';
    case {'.nii','.NII'}
        ImReader = 'nii_reader';
    case {'.vtk','.VTK'}
        ImReader = 'ReadVTK2';
    otherwise
        ImReader = 'ReadData3D';
end

bEndingSpecified = 0;
if nargin>=2
    if isfield(options,'bl'), bl = options.bl; bEndingSpecified = 1; end
    if isfield(options,'bViewIsosurface'), bViewIsosurface = options.bViewIsosurface; end
end
if exist(file,'file');
    switch ImReader
        case 'PhilipsImage'
            fprintf('Reading raw image %s (assuming Philips format) ...\n',file)
            [PH,hPH] = io_ConvertPhilipsImage(file);
            hPH = ParseHeader(hPH);
        case 'nii_reader'
            fprintf('Reading image %s (nifti format)... ',file);
            niiImage = load_untouch_nii(file);
            fprintf(' Finished!\n ');
            PH  = niiImage.img;
            hPH = ParseHeader(niiImage); 
        case 'ReadData3D'
            fprintf('Reading image %s (with ReadData3D)... ',file);
            [PH,info]=ReadData3D(file);
            fprintf(' Finished!\n ');
            if numel(find(PH))==0
                % A hack over a strange bug in ReadData3D with GIPL:
                PH = gipl_read_volume(info);
            end
            hPH = ParseHeader(info); 
        case 'ReadVTK2'
            fprintf('Reading image %s (vtk format)... ',file);
            if(bEndingSpecified)
                fprintf(' ending: %s',bl);
                [PH,hPH] = read_image_vtk2(file,bl);
            else
                [PH,hPH] = read_image_vtk2(file);
            end
            % Hack to recover unexpected behaviour in one case:
            if max(PH(:)) == 0
                if numel(unique(PH(:)))==2
                    PH = -PH;
    %                 Min = min(unique(PH));
    %                 PH = PH - Min;
                end
            end
            fprintf(' Finished!\n ');
            hPH = ParseHeader(hPH);
    end
else
    fprintf('ERROR! image file does not exist: %s\n',file);
    PH = NaN;
    hPH = NaN;
end
if(bViewIsosurface)
    show_segment_surface(PH,hPH.Mv2w);
end