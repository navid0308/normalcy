%**************************************************************************
%
%  PROGRAM TITLE      Calc_braincase.m
%
%  WRITTEN BY        Kirk E. Smith and Gregory G. Reiker
%  DATE WRITTEN      July 11, 2011 
%  WRITTEN FOR       Pediatric head modeling project
%
%  REVISIONS BY      Gregory G. Reiker, September 20, 2011
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [ bc_vol, sk_l, sk_w, sk_h, sk_geo_mean, sk_perim ] = Calc_braincase ...
%           ( g_3d, num_row, num_col, num_pln, kern, s_3d, local_z, ave_bi_edge_hu, ...
%             outputs_directory, out_file_name_root ) ;
%
%     where
%       bc_vol          braincase image volume.
%       sk_l            skull length.
%       sk_w            skull width.
%       sk_h            skull height.
%       sk_geo_mean     skull geometric mean.
%       sk_perim        skull perimeter.
%       g_3d            grayscale image volume.
%       num_row         number of rows of volume.
%       num_col         number of columns of volume.
%       num_pln         number of planes of volume.
%       kern            kernel for erosion and dilation.
%       s_3d            segmented image volume.
%       local_z         clip plane value.
%       ave_bi_edge_hu  average HU value of inside bone edge for
%                       thresholding. 
%       outputs_directory the name of the directory (folder) where
%                             the outputs will be located.
%       out_file_name_root the root name of the output files that are
%                             generated.
%
%  PROGRAM DESCRIPTION
%      Calc_braincase is a function that creates virtual endocasts from CT data
%   CT data is assumed to be 0.5 mm isotropic voxels. A series of erosions
%   and dilations are used to isolate the braincase in order to measure
%   it's volume. The threshold currently is determined based on the find
%   edges function.  The value of braincase volume is returned to
%   normalcy and written to the csv file. This function should be called
%   after the main loop of normalcy. 
%       Added in the calculation of the skull lxwxh and geometric mean.
%   local_z comes from calc_landmarks which was edited to reflect the
%   coordinates after they have been auto aligned to the nas, par, pal
%   plane.     
%
%  FILES
%       Outputs (located in "outputs_directory"):
%        .bc.vol    file with output braincase volume
%
%  DEPENDENCIES
%         MATLAB     (win64)                Version 7.12.0.635 (R2011a)
%         Image Processing Toolbox          Version 7.2        (R2011a)
%         Signal Processing Toolbox         Version 6.15       (R2011a)
%         Statistics Toolbox                Version 7.5    
%
%     Code dependencies are indicated by the level of indent:
%       calc_circumf_angles
%       rotate_ct_circumf
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------   ----------------------------------------
%       1.0     July 11, 2011     Initial release.
%
%  COPYRIGHT
%
% Copyright (c) 2011 Washington University in St. Louis
%  
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%**************************************************************************

function [ bc_vol, sk_l, sk_w, sk_h, sk_geo_mean, sk_perim ] = Calc_braincase ...
    ( g_3d, num_row, num_col, num_pln, kern, s_3d, local_z, ave_bi_edge_hu, ...
      outputs_directory, out_file_name_root )

    % Braincase segmentation
    disp('Braincase segmentation ...') ;
    skull_3d_bin = zeros(size(g_3d)) ;
%    skull_thres=225;
    skull_thres=ave_bi_edge_hu ; 
    for pln_num = 1:num_pln
        for col_num = 1:num_col
            for row_num = 1:num_row
                if g_3d(row_num, col_num, pln_num) >= skull_thres
                     skull_3d_bin(row_num, col_num, pln_num) = 1 ;
                else
                     skull_3d_bin(row_num, col_num, pln_num) = 0 ;
                end
            end
        end
    end
% remove noise
    P = 100 ; % defines size of connected pixels to remove if smaller than
    skull_3d_bin = bwareaopen(skull_3d_bin, P, 6) ;
    skull_3d_compl=imcomplement(skull_3d_bin);
    % Set the first two and last two rows, columns, and planes of the
    % segmentation image to zero. Some of the algorithms used below require
    % 2 "guard planes" on all sides. Clip at PAL-PAR plane, hard coded now.

    disp( ...
        'Clipping slices below the PAR-PAL plane ...') ;
%     skull_3d_compl(1:2,       :,         :        ) = 0 ;
%     skull_3d_compl(end-1:end, :,         :        ) = 0 ;
%     skull_3d_compl(:,         1:2,       :        ) = 0 ;
%     skull_3d_compl(:,         end-1:end, :        ) = 0 ;
    skull_3d_compl(:,         :,         1:local_z-1      ) = 0 ;
%     skull_3d_compl(:,         :,         end-1:end) = 0 ;

    disp( ...
        'Begin iterative erosion to isolate braincase ...') ;    
%    skull_3d_erode = skull_3d_compl;
    kern = logical(kern) ;
    skull_3d_erode = logical(skull_3d_compl .* s_3d);
    for n=1:10
        skull_3d_erode = imerode(skull_3d_erode, kern);
    end 
    braincase_cc = bwconncomp(skull_3d_erode, 6);
    braincase_cc_lb = labelmatrix(braincase_cc);
    c=regionprops(s_3d,'Centroid');
    cx=round(c.Centroid(1));
    cy=round(c.Centroid(2));
    cz=round(c.Centroid(3));
    brain_lb=braincase_cc_lb(cx,cy,cz);
    % It is possible that the erosions leave a hole and therefore the
    % center seed point isn't part of the object. Should modify by finding
    % nearest object to that seed point. As a temporary workaround, setting
    % the object to 1. This may or may not return the brain object, but it
    % should keep the program from crashing. Now trying biggest instead.
    if brain_lb==0
            disp( ...
        'Brain object not found at center as expected. Setting object to 1 by default ...') ;
        bc_numPixels = cellfun(@numel,braincase_cc.PixelIdxList);
        [big_obj,bo] = max(bc_numPixels);
        brain_lb=bo;
    end
    cc=braincase_cc;
    cc.NumObjects=1;
    cc.PixelIdxList=braincase_cc.PixelIdxList(brain_lb);
    brain_3d=double(labelmatrix(cc));
    disp( ...
        'Begin iterative conditional dilation to represent braincase ...') ; 
    for n=1:11
        brain_3d = imdilate(brain_3d, kern);
        brain_3d = brain_3d.*skull_3d_compl;
    end
    vox_mm=0.5;
    bc_vol=sum(sum(sum(brain_3d)))*vox_mm*vox_mm*vox_mm;
    mid_sag=rot90(squeeze(brain_3d(:,cy,:)));

    figure(9);
    imagesc(mid_sag);
    colormap('gray') ;
    title(['Sagittal plane ' num2str(cy)]) ;
    mid_cor=rot90(squeeze(brain_3d(cx,:,:)));

    figure(10);
    imagesc(mid_cor);
    colormap('gray') ;
    title(['Coronal plane ' num2str(cx)]) ;
    
    full_bcmat_file_name = [outputs_directory '\' out_file_name_root ...
        '.bc.mat'] ;
    save(full_bcmat_file_name , 'mid_sag', 'mid_cor', 'brain_3d');
    
    full_bc_file_name = [outputs_directory '\' out_file_name_root ...
        '.bc.vol'] ;
    fid=fopen(full_bc_file_name,'w');
    fwrite(fid,brain_3d,'uint8');
    fclose(fid);    

     % Calculate skull bounding box for geometric mean calculation.
     
    disp('Calculating skull lxwxh and geometric mean ...') ;
    disp( ...
        'Clipping slices below the PAL-PAR plane ...') ;
    skull_3d_bin(:,         :,         1:local_z-1      ) = 0 ;
    
    % Connect skull to remove any background. Assumes skull is largest
    % remaining object.
    disp( ...
        'Connecting skull to remove any background ...') ;
    skull_cc = bwconncomp(skull_3d_bin, 6);
        numPixels = cellfun(@numel,skull_cc.PixelIdxList);
        [biggest,idx] = max(numPixels);
        skull_3d_obj=skull_3d_bin ;
        skull_3d_obj(1:end)= 0 ;
        skull_3d_obj(skull_cc.PixelIdxList{idx}) = 1;
    
    [temp1, temp2, temp3] = ind2sub(size(skull_3d_obj), find(skull_3d_obj > 0)) ;
    list_bin_rcp = [temp1 temp2 temp3] ;
    min_rcp_of_data = [min(list_bin_rcp(:,1)) ...
                       min(list_bin_rcp(:,2)) ...
                       min(list_bin_rcp(:,3))] ;
    max_rcp_of_data = [max(list_bin_rcp(:,1)) ...
                       max(list_bin_rcp(:,2)) ...
                       max(list_bin_rcp(:,3))] ;
   
   % Skull Geometric Mean = cubed root of length x width x height, where length
   % width and height are the maximum dimensions of the bounding box in the
   % x, y, z (rcp) axes. Length = change in row, width = change in column,
   % height = change in plane
   sk_l = (max_rcp_of_data(1) - min_rcp_of_data(1)) * vox_mm ;
   sk_w = (max_rcp_of_data(2) - min_rcp_of_data(2)) * vox_mm ;
   sk_h = (max_rcp_of_data(3) - min_rcp_of_data(3)) * vox_mm ;
   sk_geo_mean = (sk_l * sk_w * sk_h)^(1/3) ;
   
    % Calculate the circumference of the skull defined by the par-pal
    % coordinate system and rotated about x axis (Analyze system) such that
    % the frontal pole (fp) and the occipital pole (op) are on the same
    % transverse plane.  min and max return the indices into the array, so
    % to get the 3d pt, find the corresponding max or min pt, but then get
    % the full 3d coordinates for that pt. 
   
    % calculate the fp and op and list in x,y,z format
    [ind_max r]= ind2sub(size(list_bin_rcp(:,1)), find(list_bin_rcp(:,1) == max_rcp_of_data(1)));
    [ind_min r]= ind2sub(size(list_bin_rcp(:,1)), find(list_bin_rcp(:,1) == min_rcp_of_data(1)));
    fp = [list_bin_rcp(ind_max(1),2) list_bin_rcp(ind_max(1),1) list_bin_rcp(ind_max(1),3)] ; 
    op = [list_bin_rcp(ind_min(1),2) list_bin_rcp(ind_min(1),1) list_bin_rcp(ind_min(1),3)] ;

    [image_center_in, alpha, beta, gamma, fp, op] = calc_circumf_angles(fp, op) ; 
        image_center_out = image_center_in ;
    disp(['Aligned CT op = ' num2str(op(3))]) ; 
    disp(['Aligned CT fp = ' num2str(fp(3))]) ;
     
    % Combine braincase segmentation with skull before rotating. This helps
    % with closing and filling in the cases of open sutures and fontanelles
    
    sk_br_obj = brain_3d + skull_3d_obj ;
    rowi = find(sk_br_obj > 0) ;
    skull_brain_bin = zeros(size(sk_br_obj)) ;
    skull_brain_bin(rowi) = 1 ;
    skull_3d_obj = skull_brain_bin ;
    % align the binary skull to the op-fp plane and return
    [skull_circumf]= rotate_ct_circumf(skull_3d_obj, image_center_in, image_center_out, ...
        alpha, beta, gamma);
   
    circumf_2d_image = skull_circumf(: , : , op(3)) ; 

    se = strel('disk',10);
    circumf_2d_image_closed = imclose(circumf_2d_image, se);
    pfill=bwfill(circumf_2d_image_closed,'holes');

    figure(11);
    imagesc(pfill);
    colormap('gray') ;
    title(['Perimeter skull fill plane ']) ;
    pc=regionprops(pfill,'Perimeter');
    sk_perim=pc.Perimeter * 0.5; %mm
    disp(['skull perimeter(mm) = ' num2str(sk_perim)]) ;
    
end

