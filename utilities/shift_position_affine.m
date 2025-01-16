function [shifted_array,indices] = shift_position_affine(X,pos1,ori1,pixdim1,sizey,pos2,ori2,pixdim2,d1,d2,int,cutoff)
    % Input X is a 3-d binary array with pixeldimensions given in pixdim1
    % The nonzero indices are the basiscoefficients of the affine transformed
    % coordinate system given by the translation in pos1 and rotation in ori1
    % 
    % Output is the same binary 3-d voxels but in the coordinate system and
    % array size given in sizey, pos2, ori2 and pixdim2
    str4='nocutoff';
    if length(size(X))==3
        str3="interpolation";
        % interpolate the size of the data in the first coordinate system onto the size of the second
        if strcmp(str3,int)
            ratio=pixdim1./pixdim2;
            x1=(size(X,1)-1)*pixdim1(1);
            x1_v=0:x1/(size(X,1)-1):x1;
            if ratio(1)>1
                x1_vq=0:x1/(round(ratio(1)*size(X,1))-1):x1;
                pixdim1(1)=x1/(round(ratio(1)*size(X,1))-1);
            else 
                x1_vq=x1_v;
            end
            x2=(size(X,2)-1)*pixdim1(2);
            x2_v=0:x2/(size(X,2)-1):x2;
            if ratio(2)>1
                x2_vq=0:x2/(round(ratio(2)*size(X,2))-1):x2;
                pixdim1(2)=x2/(round(ratio(2)*size(X,2))-1);
            else
                x2_vq=x2_v;
            end
            x3=(size(X,3)-1)*pixdim1(3);
            x3_v=0:x3/(size(X,3)-1):x3;
            if ratio(3)>1
                x3_vq=0:x3/(round(ratio(3)*size(X,3))-1):x3;
                pixdim1(3)=x3/(round(ratio(3)*size(X,3))-1);
            else
                x3_vq=x3_v;
            end
            [rect1x,rect1y,rect1z]=ndgrid(x1_v,x2_v,x3_v);
            [rect_interpx,rect_interpy,rect_interpz]=ndgrid(x1_vq,x2_vq,x3_vq);
            X_interp=interpn(rect1x,rect1y,rect1z,X,rect_interpx,rect_interpy,rect_interpz);
        else
            X_interp=X;
        end
        [indices1_x,indices1_y,indices1_z]=ind2sub(size(X_interp),find(X_interp));
        indices1=[indices1_x indices1_y indices1_z];
    elseif size(X,2)==3
        indices1=X;
    else
        error("wrong dimensions")
    end
    %store the orientation vectors in rotation matrices. The direction of z
    %axis needs to be specified
    ori1_y=[ori1(1); ori1(2); ori1(3)];
    ori1_x=[ori1(4); ori1(5); ori1(6)];
    ori1_z=d1*cross(ori1_y,ori1_x);
    ori2_y=[ori2(1); ori2(2); ori2(3)];
    ori2_x=[ori2(4); ori2(5); ori2(6)];
    ori2_z=d2*cross(ori2_y,ori2_x);
    Rot1=[ori1_x ori1_y ori1_z];
    Rot2=[ori2_x ori2_y ori2_z];
    
    indices2=(Rot2^(-1))*(Rot1*((indices1-1)'.*pixdim1)-(pos2-pos1));
    indices2=round(indices2./(pixdim2))'+1;
    shifted_array=zeros(sizey);
    k=1;
    if strcmp(str4,cutoff)
        indices=indices2;
    else
        indices=[];
        for i=1:size(indices2,1)
            if all(indices2(i,:)>[0 0 0]) && all(indices2(i,:)<[sizey(1) sizey(2) sizey(3)])
                indices(end+1,:)=indices2(i,:);
                shifted_array(indices2(i,1),indices2(i,2),indices2(i,3))=1;
                k=k+1;
            end
        end
    end
end