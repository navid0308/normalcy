%**************************************************************************
%
%  PROGRAM TITLE      seg_thresh.m
%
%  WRITTEN BY         Kirk E. Smith
%  DATE WRITTEN       November 3, 2011
%  WRITTEN FOR        Pediatric head modeling project     
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [ obj_s_3d ] = seg_thresh( obj_3d, head_thres, clip_plane ) ;
%
%     where
%       obj_s_3d    segmented volume (binary)
%       obj_3d      input volume
%       head_thres  head threshold for segmentation
%       clip_plane  clip z-plane, not currently used
%
%  PROGRAM DESCRIPTION
%       This function segments the input volume based on the threshold
%       value and returns the corresponding binary volume.
%
%   Based on code from rotate_ct.m         
%
%  FILES
%     standard input - not used
%     standard output - not used
%
%  DEPENDENCIES
%         MATLAB     (win64)                Version 7.12.0.635 (R2011a)
%         Image Processing Toolbox          Version 7.2        (R2011a)
%         Signal Processing Toolbox         Version 6.15       (R2011a)
%         Statistics Toolbox                Version 7.5    
%
%     Code dependencies are indicated by the level of indent:
%
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     November 3, 2011   Initial release.
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
%************************************************
function [ obj_s_3d ] = seg_thresh( obj_3d, head_thres, clip_plane )

 obj_3d_orig = obj_3d ;
 
 % Threshold Head
 temp0 = size(obj_3d) ;
    num_row = temp0(1) ;
    num_col = temp0(2) ;
    num_pln = temp0(3) ;
    disp('making the obj head binary ...') ;
    head_3d_bin1 = zeros(size(obj_3d)) ;
%     head_thres=100;
    % head_thres=1; % set it to one for testing a binary image
    parfor pln_num = 1:num_pln
        for col_num = 1:num_col
            for row_num = 1:num_row                 
                if obj_3d(row_num, col_num, pln_num) >= head_thres
                     head_3d_bin1(row_num, col_num, pln_num) = 1 ;
                else
                     head_3d_bin1(row_num, col_num, pln_num) = 0 ;
                end
            end
        end
    end

% Once aligned need to figure out how to segment using code below
% Connect head to remove any background. Assumes head is largest
% remaining object.
    disp( ...
        'Connecting head to remove any background ...') ;
    head_cc = bwconncomp(head_3d_bin1, 6);
        numPixels = cellfun(@numel,head_cc.PixelIdxList);
        [biggest,idx] = max(numPixels);
        head_3d_obj=head_3d_bin1 ;
        head_3d_obj(1:end)= 0 ;
        head_3d_obj(head_cc.PixelIdxList{idx}) = 1;    

 % Zero the slices in the binary image below the clip plane PAR(3)
 %  head_3d_obj(:,         :,         1:clip_plane-1      ) = 0 ;
 
 % automated and aligned grayscale (g_3d) and binary files (obj_s_3d)
    obj_s_3d = head_3d_obj ;
    
% For testing -------------------------------------
%    Display the rotated head in 3D. ----------------------------
% Something is not right with the display. x and y axes are interchanged.
% So even though data is in x,y,z MATLAB thinks x is y and y is x.
% head_3d_obj_xyz = permute(head_3d_obj,[2 1 3]);
% figure(2) ;
% p = patch(isosurface(head_3d_obj_xyz, 0.5)) ;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none') ;
% daspect([1 1 1]) ;
% view(90,0) ;
% camlight ;
% lighting gouraud ;
%------------------------------------------------------------------------


end

