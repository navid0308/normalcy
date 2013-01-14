%**************************************************************************
%
%  PROGRAM TITLE      rotate_ct_circumf.m
%
%  WRITTEN BY         Kirk E. Smith
%  DATE WRITTEN       September 14, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%        function [skull_circumf]= rotate_ct_circumf...
%           ( skull_3d_obj, image_center_in, image_center_out, alpha, beta, gamma) ;
%
%     where
%       skull_circumf      aligned skull volume (binary)
%       skull_3d_obj       skull volume (binary)
%       image_center_in    row, column, plane of input image center.
%       image_center_out   rcp of output image center.
%       alpha              angle of rotation about x, degrees
%       beta               angle of rotation about y, degrees
%       gamma              angle of rotation about z, degrees
%
%  PROGRAM DESCRIPTION
%           This function rotates the skull volume with the input center
%       and rotation angles previously calculated by calc_circumf_angles.m .  
%
%       This function is based on rotate_3d_ct.m
%       (based on http://blogs.mathworks.com/steve/2006/08/17/spatial-
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
%       1.0     September 14, 2011 Initial release.
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
function [skull_circumf]= rotate_ct_circumf...
    ( skull_3d_obj, image_center_in, image_center_out, alpha, beta, gamma)

% Make a 3D affine tform structure.
% T0 --> subtract out the center of rotation in the input space
% T1 --> rotate by Euler angle alpha
% T2 --> rotate by Euler angle beta
% T3 --> rotate by Euler angle gamma
% T4 --> add in the center of rotation in the output space

% % For testing -------------------------------------
% % Display original Binary image
% 
% figure ;
% p = patch(isosurface(head_3d_bin, 0.5)) ;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none') ;
% daspect([1 1 1]) ;
% view(90,0) ;
% camlight ;
% lighting gouraud ;
% %---------------------------------------------------------------

% Set g_3d to the binary skull so code doesn't have to be modified.

 g_3d = skull_3d_obj ;

% interchanging g_3d x and y analyze coordinates to MATLAB rcp.
g_3d = permute(g_3d,[2 1 3]);

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

tformfwd(image_center_in, tform) ;

% Create all the input arguments needed for tformarray, which is the Matlab
% function that does the work of rotating the image.

R = makeresampler('linear', 'fill') ;
TDIMS_A = [1 2 3] ;
TDIMS_B = [1 2 3] ;
TSIZE_B = size(g_3d) ;
TMAP_B  = [] ;
F = 0 ; % This value is the fill value and it is set to 0 since it is a
% binary object.

% Rotate the 3D image.
g_3d = tformarray(g_3d, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F) ;

% Interchanging g_3d back to x and y analyze coordinates from MATLAB rcp,
% since subsequent code is expecting this
g_3d = permute(g_3d,[2 1 3]);
    
% Once aligned need to figure out how to segment using code below
% Connect head to remove any background. Assumes head is largest
% remaining object.
%     disp( ...
%         'Connecting head to remove any background ...') ;
%     head_cc = bwconncomp(head_3d_bin1, 6);
%         numPixels = cellfun(@numel,head_cc.PixelIdxList);
%         [biggest,idx] = max(numPixels);
%         head_3d_obj=head_3d_bin1 ;
%         head_3d_obj(1:end)= 0 ;
%         head_3d_obj(head_cc.PixelIdxList{idx}) = 1;    

 % Zero the slices in the binary image below the clip plane PAR(3)
%  head_3d_obj(:,         :,         1:image_center_in(3)-1      ) = 0 ;
 
% %    Display the rotated head in 3D. ------------------------------------
% % head_3d_obj = head_3d_bin1; % only used if connected components commented out
% figure ;
% p = patch(isosurface(head_3d_obj, 0.5)) ;
% set(p, 'FaceColor', 'red', 'EdgeColor', 'none') ;
% daspect([1 1 1]) ;
% view(90,0) ;
% camlight ;
% lighting gouraud ;
% %------------------------------------------------------------------------

% automated and aligned grayscale (g_3d) and binary files (s_3d)
skull_circumf = g_3d ;

end