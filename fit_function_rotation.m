%**************************************************************************
%
%  PROGRAM TITLE      fit_function_rotation.m
%
%  WRITTEN BY         Gregory G. Reiker
%  DATE WRITTEN       September 29, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY       Kirk Smith 
%  DATE MODIFIED      November 1, 2011
%
%  CALLING SYNTAX
%     Use the following syntax:
%        function [sse] = ...
%           fit_function_rotation(parms_in, ref_point_mat_xyz, mr_list_border_bc_rcp, ...
%                       image_center_in, local_origin, step_size) ;
%       in code,
%           f = @(x)fit_function_rotation(x, cum_ct_pts_xyz, ...
%                        mr_list_border_bc_rcp, center_xyz, ...
%                        local_origin, step_size) ;
%           [parms_out, sse_reg, exitflag, foutput] = fminsearch(f, parms_init, optimset('MaxFunEvals',4000)) ;
%
%     where
%       sse                     sum of squared errors between estimated MR and aligned CT input
%       parms_in                input parameter:  image center out vector,
%                                                  alpha, beta, gamma rotation angles
%       ref_point_mat_xyz       reference points matrix (xyz)
%       mr_list_border_bc_rcp   MR list of bounding box border before cleanup (rcp)
%       image_center_in         input image center for rotation
%       local_origin            volume local origin
%       step_size               step size for iterations
%
%  PROGRAM DESCRIPTION
%       Estimate the parameters of the best-fit rotation and translation
%     using fminsearch, passing in MR border points. 
%     Calculate corresponding points within the function.       
%
%  FILES
%     standard input - not used
%     standard output - not used
%
%  DEPENDENCIES
%         MATLAB     (win64)                Version 7.12.0.635 (R2011a)
%         Image Processing Toolbox          Version 7.2        (R2011a)
%         Optimization Toolbox              Version 6.0        (R2011a)
%         Signal Processing Toolbox         Version 6.15       (R2011a)
%         Statistics Toolbox                Version 7.5    
%
%     Code dependencies are indicated by the level of indent:
%
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     September 29, 2011 Initial release.
%       1.1     November 1, 2011   New approach: Minimum distance from
%                                  reference point to boundary.
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

% Pass in MR boundary points
function [sse] = ...
    fit_function_rotation(parms_in, ref_point_mat_xyz, mr_list_border_bc_rcp, ...
                          image_center_in, local_origin, step_size)
    
    image_center_out(1) = parms_in(4) + image_center_in(1)-step_size ;
    image_center_out(2) = parms_in(5) + image_center_in(2)-step_size ;
    image_center_out(3) = parms_in(6) + image_center_in(3)-step_size ;
    alpha     = parms_in(1)-step_size ;
    beta     = parms_in(2)-step_size ;
    gamma     = parms_in(3)-step_size ;
    
   % translation
T0 = [1     0     0     0 ; ...
      0     1     0     0 ; ...
      0     0     1     0 ; ...
      -image_center_in  1] ;
  
% T1, T2, T3 are about x, y, z in a left handed coordinate system

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

    %Apply rotations and translations
    tform = maketform('affine', T) ;
    %change from rcp to xyz
    mr_list_border_bc_xyz = mr_list_border_bc_rcp(:, [2 1 3]);
%    mr_est_xyz = tformfwd(mr_ref_point_mat_xyz, tform) ;
    mr_list_border_bc_xyz = tformfwd(mr_list_border_bc_xyz, tform) ;
%    mr_par = tformfwd(mr_par, tform) ;
%
%_______________________________________________________________________
% This code calculates the minimum distance from the ct landmark to the mr
% surface boundary.
border_x = mr_list_border_bc_xyz( : ,1) ;
border_y = mr_list_border_bc_xyz( : ,2) ;
border_z = mr_list_border_bc_xyz( : ,3) ;
landmark_x = ref_point_mat_xyz( : ,1) ;
landmark_y = ref_point_mat_xyz( : ,2) ;
landmark_z = ref_point_mat_xyz( : ,3) ;

x1 =  local_origin(1);
y1 = local_origin(2) ;
z1 = local_origin(3) ;
x2 = landmark_x ;
y2 = landmark_y ;
z2 = landmark_z ;
local_z = z1;

% calulate the lowest z value in ct ref points and set the clip half way
% between that and local origin
%[clip_v, clip_i]  = min(landmark_z) ; 
%clip_z = (clip_v - local_z)/2 ;
clip_z = 1;
% Find the border voxels excluding the ones less than local_z
[index] = find(border_z > clip_z) ;
border_x = border_x(index) ;
border_y = border_y(index) ;
border_z = border_z(index) ;
mr_list_border_bc_xyz_thin = [border_x border_y border_z] ;

temp1 = size(ref_point_mat_xyz) ;
land_tot = temp1(1) ;
mr_point_mat_xyz = zeros(land_tot, 3) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New approach. Minimum distance from ref pt to boundary

for land = 1:land_tot
%disp('Finding minimum distance ref pt to boundary ...') ;
        dist_ctref_mrbound = ...
            ((mr_list_border_bc_xyz_thin(:,1) - ref_point_mat_xyz(land, 1)).^2) ...
               + ((mr_list_border_bc_xyz_thin(:,2) - ref_point_mat_xyz(land, 2)).^2) ...
               + ((mr_list_border_bc_xyz_thin(:,3) - ref_point_mat_xyz(land, 3)).^2) ;
        [~, short_index] = min(dist_ctref_mrbound) ;
  
    mr_point_mat_xyz(land, : ) = [border_x(short_index)...
        border_y(short_index)...
        border_z(short_index)] ;
    
%     mr_point_mat_xyz(land, : ) = [border_x(spos_index(short_index))...
%         border_y(spos_index(short_index))...
%         border_z(spos_index(short_index))] ;
    
end
temp = size(ref_point_mat_xyz) ;
    num_ref_points = temp(1) ;
    if (num_ref_points <= 0)
        disp('ERROR:  NO REFERENCE POINTS SPECIFIED. EXITING!') ;
        return ;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%_______________________________________________________________________

    % sum of errors squared between estimated MR and aligned CT input
    sse = sum(((mr_point_mat_xyz(:,1) - ref_point_mat_xyz(:,1)).^2) ...
           + ((mr_point_mat_xyz(:,2) - ref_point_mat_xyz(:,2)).^2) ...
           + ((mr_point_mat_xyz(:,3) - ref_point_mat_xyz(:,3)).^2)) ;

 end