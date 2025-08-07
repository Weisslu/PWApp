 % Author: LWeissinger, HVoss, 20.4.2023-18.07.2025
% Read dicom directory with phase contrast data
% The dicom data need to be downloaded from Ambra (not xnat or nile), only then the file order makes sense
function [tdim,res,pixdim,delays,phis,mags,orientation,position,directionFlag]=read_dicoms()
f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
[~,directory] = uigetfile('*.dcm','Select dicom data');
delete(f); %delete the dummy figure
if directory ~= 0
    files = dir([directory '/*.dcm']);
    test_num = str2double(directory((end-2):(end-1)));
    numberoffiles = size(files,1);
    %slice_positions = zeros(numberoffiles,1);
    filename = fullfile(directory, files(1).name);
    info = dicominfo(filename);
    tdim = info.CardiacNumberOfImages;
    xdim = info.Rows;  
    ydim = info.Columns;
    zdim = 0;
    innerdim = 'slice';
    fastline = 0;
    i = 1; k = 1;
    while i <= numberoffiles
        filename = fullfile(directory, files(i).name);
        info = dicominfo(filename);
        slice_positions(k) = info.SliceLocation;
        zdim_old = zdim;    
        zdim = length(unique(slice_positions));
        if zdim_old == zdim
            if i == 2
                innerdim = 'time';
                fastline = 1;
            end
            if strcmp(innerdim,'slice') 
                i = numberoffiles;
            end
        end
        if fastline == 1
            i = i + tdim - 1;
        end
        i = i + 1;
        k = k + 1; %ensures slice_positions is not filled with zeros
    end
    filename = fullfile(directory, files(1).name);
    info = dicominfo(filename);
    %zdim=info.Private_0021_104f;%lol.
    res = [xdim; ydim; zdim];
    if exist('info.SpacingBetweenSlices','var') % Implement via slice_positions
        pix_z = info.SpacingBetweenSlices;
    else
        pix_z = info.SliceThickness;
    end
    pixdim=[info.PixelSpacing ; pix_z];
    orientation = info.ImageOrientationPatient;
    position = info.ImagePositionPatient;
    
    filename_phis = sprintf('phis_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
    filename_mags = sprintf('mags_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
    filename_delays = sprintf('delays_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
    h = waitbar(0,"Loading velocity data");
    set(h,'Pointer','watch');
    drawnow()

    mags = zeros(xdim,ydim,zdim,tdim);
    phis = mags;
    phis_x = mags;
    phis_y = mags;
    phis_z = mags;
    delays = zeros(zdim,tdim);
    if numberoffiles ~= 2*tdim*zdim && numberoffiles ~= 5*tdim*zdim && numberoffiles ~= 4*tdim*zdim
        error('Something is wrong, probably not velocity data that you tried to load!') 
        %directory=0;tdim=0;res=0;pixdim=0;delays=0;phis=0;mags=0;vangio3d=0;orientation=0;position=0;directionFlag=0;
    elseif numberoffiles == 2*tdim*zdim || numberoffiles == 5*tdim*zdim
        done = 0;
        fileind = 0;
        for twice = 1:2
            for z = 1:zdim
                for t = 1:tdim
                    if strcmp(innerdim,'slice')
                        fileind = (t-1)*(zdim) + (z) + (twice-1)*tdim*zdim;
                    elseif strcmp(innerdim,'time')
                	    fileind = fileind+1;
                    end
                    if isvalid(h)
                        waitbar(((z-1)*(tdim) + (t) + (twice-1)*tdim*zdim)/double(2*tdim*zdim),h,"Loading velocity data");
                    end
                    % Images
            
                    filename = fullfile(directory, files(fileind).name);
                    data = dicomread(filename);
                    if z == 2 && done == 0
                        info2 = dicominfo(filename);
                        position2 = info2.ImagePositionPatient;
                        ori_x = [orientation(1); orientation(2); orientation(3)];
                        ori_y = [orientation(4); orientation(5); orientation(6)];
                        ori_z = cross(ori_y,ori_x);
                        pos2 = position + pixdim(3)*ori_z;
                        ori_z = -ori_z;
                        pos3 = position + pixdim(3)*ori_z;
                        if norm(position2-pos2) < 0.2
                            directionFlag = -1;
                        elseif norm(position2-pos3) < 0.2
                            directionFlag = 1;
                        else
                            error("most likely there is an implementation error")
                        end
                        done = 1;
                    end
                    %  image(data); colormap gray; title(num2str(fileind));
            
                    if fileind <= tdim*zdim
                        phis(:,:,z,t) = data;
                        info = dicominfo(filename);                   
                        delays(z,t) = info.TriggerTime;
                    else
                        mags(:,:,z,t) = data;
                    end 
                end
            end
        end
    elseif numberoffiles == 4*tdim*zdim
        done = 0;
        fileind = 0;
        for twice = 1:4
            for z = 1:zdim
                for t = 1:tdim
                    if strcmp(innerdim,'slice')
                        fileind = (t-1)*(zdim) + (z) + (twice-1)*tdim*zdim;
                    elseif strcmp(innerdim,'time')
                	    fileind = fileind+1;
                    end
                    if isvalid(h)
                        waitbar(((z-1)*(tdim) + (t) + (twice-1)*tdim*zdim)/double(4*tdim*zdim),h,"Loading velocity data");
                    end
                    % Images
            
                    filename = fullfile(directory, files(fileind).name);
                    data = dicomread(filename);
                    if z==2 && done==0
                        info2 = dicominfo(filename);
                        position2 = info2.ImagePositionPatient;
                        ori_x = [orientation(1); orientation(2); orientation(3)];
                        ori_y = [orientation(4); orientation(5); orientation(6)];
                        ori_z = cross(ori_y,ori_x);
                        pos2 = position+pixdim(3)*ori_z;
                        ori_z = -ori_z;
                        pos3 = position+pixdim(3)*ori_z;
                        if norm(position2-pos2) < 0.2
                            directionFlag = -1;
                        elseif norm(position2-pos3) < 0.2
                            directionFlag = 1;
                        else
                            error("most likely there is an implementation error")
                        end
                        done = 1;
                    end
                    %  image(data); colormap gray; title(num2str(fileind));
            
                    if fileind <= tdim*zdim
                        phis_x(:,:,z,t) = data;
                        info = dicominfo(filename);                   
                        delays(z,t) = info.TriggerTime;
                    elseif fileind <= 2*tdim*zdim
                        phis_y(:,:,z,t) = data;
                    elseif fileind <= 3*tdim*zdim
                        phis_z(:,:,z,t) = data;
                    else
                        mags(:,:,z,t) = data;
                    end 
                end
            end
        end
        phis = (phis_x.^2+phis_y.^2+phis_z.^2).^(1/2);
    end
delete(h);
else
    directory = 0; tdim = 0; res = 0; pixdim = 0; delays = 0; phis = 0; mags = 0; orientation = 0; position = 0; directionFlag = 0;
end
end


