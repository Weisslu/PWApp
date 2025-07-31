function [directory,tdim,res,pixdim,delays,phis,mags,vangio3d,orientation,position,directionFlag] = ...
    loadPARREC()
% Get and load input directory
f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]); %create a dummy figure so that uigetfile doesn't minimize our GUI
[filename,directory] = uigetfile('*.rec','Select .rec Data in each axis', 'MultiSelect', 'on');
%fBase = filename(1:end-5);
delete(f); %delete the dummy figure
warning('off','all');
%% grab each parrec and save corresponding data
%disp('Loading data')
PARRECFILE = fullfile(directory,filename{end-2});
[IMG1,~] = readrec_V4_2(PARRECFILE, 'noscale');
IMG1 = double(IMG1);
if size(IMG1,6)>1
    vx = squeeze(IMG1(:,:,:,:,:,2,:))-2048;
    mag1 = squeeze(IMG1(:,:,:,:,:,1,:));
else
    vx = squeeze(IMG1(:,:,:,:,:,1,:))-2048;
end

PARRECFILE = fullfile(directory,filename{end-1});
[IMG2,~] = readrec_V4_2(PARRECFILE,'noscale');
IMG2 = double(IMG2);
if size(IMG2,6)>1
    vy = squeeze(IMG2(:,:,:,:,:,2,:))-2048;
    mag2 = squeeze(IMG2(:,:,:,:,:,1,:));
else
    vy = squeeze(IMG2(:,:,:,:,:,1,:))-2048;
end

PARRECFILE = fullfile(directory,filename{end});
[IMG3,header] = readrec_V4_2(PARRECFILE, 'noscale');
IMG3 = double(IMG3);
if size(IMG3,6)>1
    vz = squeeze(IMG3(:,:,:,:,:,2,:))-2048;
    mag3 = squeeze(IMG3(:,:,:,:,:,1,:));
else
    vz = squeeze(IMG3(:,:,:,:,:,1,:))-2048;
end
warning('on','all');

if length(filename) == 4
    PARRECFILE = fullfile(directory,filename{1});
    [IMG4,header] = readrec_V4_2(PARRECFILE, 'noscale');
    IMG4 = double(IMG4);
    mags = squeeze(IMG4(:,:,:,:,:,1,:));
elseif length(filename) == 3
    mags = mean(cat(5,mag1,mag2,mag3),5);
else
    error("Something does not add up here")
end

vy = -vy;       % flip vy!
v = cat(5,vx,vy,vz); v = permute(v, [1 2 3 5 4]);
clear mag1 mag2 mag3 IMG1 IMG2 IMG3 vx vy vz

tdim = header.nphases;                                       % number of reconstructed frames
delays_temp = header.tbl(:,header.tblcols.ttime);
if length(filename) == 3
    delays_temp=reshape(delays_temp,[],tdim,2);
    delays=delays_temp(:,:,1);% temporal resolution, in ms
else
    delays_temp=reshape(delays_temp,2,[],tdim);
    delays=squeeze(delays_temp(1,:,:));
end
fov = header.fov;                                               % Field of view in cm
res = round([header.nrows header.ncols header.nslices]);        % number of pixels in row,col,slices
VENC = max(header.pevelocity)*10;                               % venc, in mm/s
pixdim = header.pixdim';                                         % the reconstructed resolution
ori = header.tbl(1,26);                                         % orientation number (1 - axial, 2 - sagittal, 3 - coronal)
position=[0 0 0];
orientation=eye(3);
directionFlag=1;
% scale velocity
v = (v./2048)*VENC;

% take the means
vMean = mean(v,5);
mags = mags./max(mags(:));


[~, vangio3d] = calc_angio(mags, v, VENC);
phis = squeeze(sqrt( sum( v.^2,4)));
disp('Data loading finished');
return