% PC1, 20.4.2023
% Read dicom directory with phase contrast data
% The dicom data need to be downloaded from Ambra (not xnat or nile), only then the file order makes sense
function [directory,tdim,res,pixdim,delays,phis,mags,vangio3d,orientation,position,directionFlag]=read_dicoms()
f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
[~,directory] = uigetfile('*.dcm','Select dicom data');
delete(f); %delete the dummy figure
files=dir([directory '/*.dcm']);
test_num = str2double(directory((end-2):(end-1)));
numberoffiles=size(files,1);
filename = fullfile(directory, files(1).name);
info = dicominfo(filename);
zdim=info.Private_0021_104f;%lol.
tdim=info.CardiacNumberOfImages;
xdim=info.Width; 
ydim=info.Height;
res=[xdim;ydim;zdim];
pixdim=[info.PixelSpacing ; info.SpacingBetweenSlices ];
orientation=info.ImageOrientationPatient;
position=info.ImagePositionPatient;

filename_phis = sprintf('phis_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
filename_mags = sprintf('mags_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
filename_delays = sprintf('delays_test%i_%ix_%iy_%iz_%it.mat',test_num,xdim,ydim,zdim,tdim);
h=waitbar(0,"Loading velocity data");
set(h,'Pointer','watch');
drawnow()
if 0%exist(filename_phis,"file") && exist(filename_mags,"file") && exist(filename_delays,"file")
    wbch = allchild(h);
    jp = wbch(1).JavaPeer;
    jp.setIndeterminate(1)
    load(filename_phis)
    load(filename_mags)
    load(filename_delays)
    fprintf("Data used from previous saving\n")
    filename = fullfile(directory, files(tdim+1).name);
    info2 = dicominfo(filename);
    position2=info2.ImagePositionPatient;
    ori_x=[orientation(1); orientation(2); orientation(3)];
    ori_y=[orientation(4); orientation(5); orientation(6)];
    ori_z=cross(ori_y,ori_x);
    pos2=position+pixdim(3)*ori_z;
    ori_z=-ori_z;
    pos3=position+pixdim(3)*ori_z;
    if norm(position2-pos2)<0.1
        directionFlag=-1;
    elseif norm(position2-pos3)<0.1
        directionFlag=1;
    else
        error("most likely there is an implementation error")
    end
else
    mags=zeros(xdim,ydim,zdim,tdim);
    phis=mags;
    delays=zeros(zdim,tdim);
    if numberoffiles ~= 2*tdim*zdim && numberoffiles ~= 5*tdim*zdim; 
        errordlg('Something is wrong, probably not velocity data that you tried to load!'); 
        directory=0;tdim=0;res=0;pixdim=0;delays=0;phis=0;mags=0;vangio3d=0;orientation=0;position=0;directionFlag=0;
    else
        done=0;
        fileind=0;
        for twice=1:2
            for z=1:zdim
                for t=1:tdim
                    fileind=fileind+1;
                    waitbar(fileind/double(2*tdim*zdim),h,"Loading velocity data");
                    % Images
            
                    filename = fullfile(directory, files(fileind).name);
                    data=dicomread(filename);
                    if z==2 && done==0
                        info2 = dicominfo(filename);
                        position2=info2.ImagePositionPatient;
                        ori_x=[orientation(1); orientation(2); orientation(3)];
                        ori_y=[orientation(4); orientation(5); orientation(6)];
                        ori_z=cross(ori_y,ori_x);
                        pos2=position+pixdim(3)*ori_z;
                        ori_z=-ori_z;
                        pos3=position+pixdim(3)*ori_z;
                        if norm(position2-pos2)<0.2
                            directionFlag=-1;
                        elseif norm(position2-pos3)<0.2
                            directionFlag=1;
                        else
                            error("most likely there is an implementation error")
                        end
                        done=1;
                    end
                    %  image(data); colormap gray; title(num2str(fileind));
            
                    if fileind<=tdim*zdim
                        phis(:,:,z,t)=data;
                        info = dicominfo(filename);                   
                        delays(z,t)=info.TriggerTime;
                    else
                        mags(:,:,z,t)=data;
                    end 
                end
            end
        end
        %save(filename_phis, 'phis','-v7.3');
        %save(filename_mags, 'mags','-v7.3');
        %save(filename_delays,'delays','-v7.3');
        mm = mean(mags,4); vMean = squeeze(mean(phis,4));
        vangio3d = mm.*sin( pi/2*rescale(vMean,-1,1));
    end
end
close(h);
end


