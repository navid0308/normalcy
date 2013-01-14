%**************************************************************************
%
%  PROGRAM TITLE      rotate_3d_ct.m
%
%  WRITTEN BY         David G. Politte and Kirk E. Smith
%  DATE WRITTEN       August 3, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY       Gregory G. Reiker
%  DATE MODIFIED      August 4-18, 2011
%  REVISIONS BY       Kirk E. Smith
%  DATE MODIFIED      August 5, 2011
%
%  CALLING SYNTAX
%     Use the following syntax:
%        function [g_3d_orig, s_3d, g_3d]= rotate_3d_ct( g_3d, image_center_in, image_center_out, ...
%                                                        alpha, beta, gamma) ;
%
%     where
%       g_3d_orig          original CT grayscale volume
%       s_3d               rotated and segmented CT volume (binary)
%       g_3d               output- rotated CT grayscale volume
%       g_3d               input - original CT grayscale volume
%       image_center_in    row, column, plane of input image center.
%       image_center_out   rcp of output image center.
%       alpha              angle of rotation about x, degrees
%       beta               angle of rotation about y, degrees
%       gamma              angle of rotation about z, degrees
%
%  PROGRAM DESCRIPTION
%           This function rotates the grayscale CT image with the input centers
%       and rotation angles previously calculated by Calc_rot_angles.m .  
%       Rotate_3d_ct.m thresholds the head to generate the aligned segmented
%       volume.  This function does a connect operation so stray voxels should be gone,
%       but it does not do a fill operation.  Also, it clips the segmented
%       volume at the z plane of the output center.
%
% (based on http://blogs.mathworks.com/steve/2006/08/17/spatial-
%                  transformations-three-dimensional-rotation/)
%       rot3d.m  - Demonstration of 3D image rotation   
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
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     August 3, 2011     Initial release.
%       1.1     August 4, 2011     Made into a function.
%       1.2     August 18, 2011    Added head segmentation and clipping.
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
function [g_3d_orig, s_3d, g_3d]= rotate_3d_ct( g_3d, image_center_in, image_center_out, ...
    alpha, beta, gamma)

% Make a 3D affine tform structure.
% T0 --> subtract out the center of rotation in the input space
% T1 --> rotate by Euler angle alpha
% T2 --> rotate by Euler angle beta
% T3 --> rotate by Euler angle gamma
% T4 --> add in the center of rotation in the output space

 g_3d_orig = g_3d ;

% interchanging g_3d MATLAB rcp to x and y analyze coordinates.
g_3d = permute(g_3d,[2 1 3]);

% % For testing -------------------------------------
% % Display original image
%  % Threshold Head
%  temp0 = size(g_3d) ;
%     num_row = temp0(1) ;
%     num_col = temp0(2) ;
%     num_pln = temp0(3) ;
%     disp('making the head binary ...') ;
%     head_3d_bin = zeros(size(g_3d)) ;
%      head_thres=-350;
%     % head_thres=1; % set it to one for testing a binary image
%     for pln_num = 1:num_pln
%         for col_num = 1:num_col
%             for row_num = 1:num_row                 
%                 if g_3d(row_num, col_num, pln_num) >= head_thres
%                      head_3d_bin(row_num, col_num, pln_num) = 1 ;
%                 else
%                      head_3d_bin(row_num, col_num, pln_num) = 0 ;
%                 end
%             end
%         end
%     end
% 
%     disp( ...
%         'Connecting head to remove any background ...') ;
%     head_cc = bwconncomp(head_3d_bin, 6);
%         numPixels = cellfun(@numel,head_cc.PixelIdxList);
%         [biggest,idx] = max(numPixels);
%         head_3d_obj_orig=head_3d_bin ;
%         head_3d_obj_orig(1:end)= 0 ;
%         head_3d_obj_orig(head_cc.PixelIdxList{idx}) = 1;    
% 
%  % Zero the slices in the binary image below the clip plane PAR(3)
%  head_3d_obj_orig(:,         :,         1:image_center_out(3)-1      ) = 0 ;
%      
%     
% figure ;
% p = patch(isosurface(head_3d_obj_orig, 0.5)) ;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none') ;
% daspect([1 1 1]) ;
% view(90,90) ;
% camlight ;
% lighting gouraud ;
% %---------------------------------------------------------------

% translation
T0 = [1     0     0     0 ; ...
      0     1     0     0 ; ...
      0     0     1     0 ; ...
      -image_center_in  1] ;
  
% Ti, T2, T3 are about x, y, z in a left handed coordinate system

% changed signs of sind's for T1 to match alpha rotation
T1 = [1 0            0           0  ; ...
      0 cosd(alpha) sind(alpha) 0  ; ...
      0 -sind(alpha) cosd(alpha) 0  ; ...
      0 0            0           1] ;

T2 = [cosd(beta) 0 -sind(beta) 0  ; ...
      0          1  0          0  ; ...
      sind(beta) 0  cosd(beta) 0  ; ...
      0          0  0          1] ;
   
T3 = [cosd(gamma) sind(gamma) 0 0  ; ...
      -sind(gamma)  cosd(gamma) 0 0  ; ...
      0            0           1 0  ; ...
      0            0           0 1] ;
   
T4 = [1     0     0     0  ; ...
      0     1     0     0  ; ...
      0     0     1     0  ; ...
      image_center_out  1] ;
  
% rotating about z, then y, then x;
T = T0 * T3 * T2 * T1 * T4 ; % Accounts for all rotations and translations.

tform = maketform('affine', T) ;

%tformfwd(image_center_in, tform) ;

% Create all the input arguments needed for tformarray, which is the Matlab
% function that does the work of rotating the image.

R = makeresampler('linear', 'fill') ;
TDIMS_A = [1 2 3] ;
TDIMS_B = [1 2 3] ;
TSIZE_B = size(g_3d) ;
TMAP_B  = [] ;
F = -1024 ; % The value to be used outside the boundaries of the input array.
% For CT data, set F = -1024. for Siemens.

% Rotate the 3D image.
g_3d = tformarray(g_3d, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F) ;

% Interchanging g_3d back to MATLAB rcp from x and y analyze coordinates,
% since subsequent code is expecting this
g_3d = permute(g_3d,[2 1 3]);

 % Threshold Head
 temp0 = size(g_3d) ;
    num_row = temp0(1) ;
    num_col = temp0(2) ;
    num_pln = temp0(3) ;
    disp('making the head binary ...') ;
    head_3d_bin1 = zeros(size(g_3d)) ;
    head_thres=-350;
    % head_thres=1; % set to 1 if already a binary
    parfor pln_num = 1:num_pln
        for col_num = 1:num_col
            for row_num = 1:num_row
                if g_3d(row_num, col_num, pln_num) >= head_thres
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
 head_3d_obj(:,         :,         1:image_center_out(3)-1      ) = 0 ;
 
%    Display the rotated head in 3D. ------------------------------------
% head_3d_obj = head_3d_bin1; % only used if connected components commented out
% Something is not right with the display. x and y axes are interchanged.
% So even though data is in x,y,z MATLAB thinks x is y and y is x.
%-----------------------------------------------------------------------
% Uncomment below to Display 3D head 
% head_3d_obj_xyz = permute(head_3d_obj,[2 1 3]);
% figure(1) ;
% p = patch(isosurface(head_3d_obj_xyz, 0.5)) ;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none') ;
% daspect([1 1 1]) ;
% view(0,0) ;
% camlight ;
% lighting gouraud ;
%------------------------------------------------------------------------

% automated and aligned grayscale (g_3d) and binary files (s_3d) in rcp
% format
s_3d = head_3d_obj ;


end