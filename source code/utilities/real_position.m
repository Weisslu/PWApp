function [shifted_array] = shift_position(x,pos1,ori1,pixdim1,pos2,ori2,pixdim2)
% Input is an array with 3-d indices of the nonzero values of a binary 3-d array,
% the real position of the first point in the array (pos) 
% and the orientation of the row and columns of the array via direction cosines (ori)
% and the pixel-dimensions
% Output is the real position of the corresponding indices
    indices=ind2sub(find(x));
    ori_x=[ori(1) ori(2) ori(3)];
    ori_y=[ori(4) ori(5) ori(6)];
    ori_z=cross(ori_x,ori_y);
    real=pos+(x-1).*ori_x*pixdim(1)+(x-1).*ori_y*pixdim(2)+(x-1).*ori_z*pixdim(3);
end