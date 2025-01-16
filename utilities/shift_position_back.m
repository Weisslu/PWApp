function [shifted_array] = shift_position_back(X,pos1,ori1,pixdim1,sizey,pos2,ori2,pixdim2)
% Input is an array with 3-d indices of the nonzero values of a binary 3-d array,
% the real position of the first point in the array (pos) 
% and the orientation of the row and columns of the array via direction cosines (ori)
% and the pixel-dimensions
% Output is the real position of the corresponding indices
        % interpolate the velocity data onto the angio-gridsize
        ratio=pixdim1./pixdim2;
        x1=(size(X,1)-1)*pixdim1(1);
        x1_v=0:x1/(size(X,1)-1):x1;
        x1_vq=0:x1/(round(ratio(1)*size(X,1))-1):x1;
        x2=(size(X,2)-1)*pixdim1(2);
        x2_v=0:x2/(size(X,2)-1):x2;
        x2_vq=0:x2/(round(ratio(2)*size(X,2))-1):x2;
        x3=(size(X,3)-1)*pixdim1(3);
        x3_v=0:x3/(size(X,3)-1):x3;
        x3_vq=0:x3/(round(ratio(3)*size(X,3))-1):x3;
        [rect1x,rect1y,rect1z]=ndgrid(x1_v,x2_v,x3_v);
        [rect_interpx,rect_interpy,rect_interpz]=ndgrid(x1_vq,x2_vq,x3_vq);
        pixdim1(1)=x1/(round(ratio(1)*size(X,1))-1);
        pixdim1(2)=x2/(round(ratio(2)*size(X,2))-1);
        pixdim1(3)=x3/(round(ratio(3)*size(X,3))-1);
        X_interp=interpn(rect1x,rect1y,rect1z,X,rect_interpx,rect_interpy,rect_interpz);
    ori1_y=[ori1(1) ori1(2) ori1(3)];
    ori1_x=[ori1(4) ori1(5) ori1(6)];
    ori1_z=cross(ori1_y,ori1_x);
    ori2_y=[ori2(1) ori2(2) ori2(3)];
    ori2_x=[ori2(4) ori2(5) ori2(6)];
    ori2_z=cross(ori2_y,ori2_x);
    [indices1_x,indices1_y,indices1_z]=ind2sub(size(X_interp),find(X_interp));
    %indices1=[indices1_x indices1_y indices1_z];
    pos_real=pos1'+(indices1_x-1).*ori1_x*pixdim1(1)+(indices1_y-1).*ori1_y*pixdim1(2)-(indices1_z-1).*ori1_z*pixdim1(3);
    
    indices2_x=round((pos_real-pos2')./(ori2_x*pixdim2(1))+1);
    indices2_y=round((pos_real-pos2')./(ori2_y*pixdim2(2))+1);
    indices2_z=round((pos_real-pos2')./(ori2_z*pixdim2(3))+1);
    shifted_array=zeros(sizey);
    indices2=[indices2_x(:,2) indices2_y(:,1) indices2_z(:,3)];
    k=1;
    for i=1:size(indices2,1)
        if all(indices2(i,:)>[0 0 0]) && all(indices2(i,:)<[sizey(1) sizey(2) sizey(3)])
            shifted_array(indices2(i,1),indices2(i,2),indices2(i,3))=1;
            k=k+1;
        end
    end
    end