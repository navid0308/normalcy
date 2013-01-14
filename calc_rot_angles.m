%**************************************************************************
%
%  PROGRAM TITLE      calc_rot_angles.m
%
%  WRITTEN BY         Gregory G. Reiker and Kirk Smith
%  DATE WRITTEN       August 4, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY      
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [image_center_in, alpha, beta, gamma, nas, par, pal ] = ...
%            calc_rot_angles( full_ref_local_file_name ) ;
%
%     where
%       image_center_in row, column, plane of image center.
%       alpha           angle of rotation about x.
%       beta            angle of rotation about y.
%       gamma           angle of rotation about z.
%       nas                 NAS
%       par                 PAR
%       pal                 PAL
%       full_ref_local_file_name name of the text file that
%                             contains the (x,y,z) coordinates of the
%                             CT reference points (NAS,PAR,PAL) that define a
%                             local coordinate system, with the name of the 
%                             directory (folder) where the inputs are located.
%
%  PROGRAM DESCRIPTION
%       Calculate a Euler angles to align the PAR, PAL, NAS coordinates
%   automatically. When applying the rotations, order of execution matters.
%   Make sure to follow it. First step is to calculate the rotations in
%   degrees. We will use the PAR as the reference.
%
%   8/10/2011 - angles seem to be calculated correctly.
%   Shift the PAR point to 0 by subtracting it from both PAL and
%   NAS and then calculating new points that way.    
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
%          read_ref_landmarks_file 
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     August 4, 2011   Initial release.
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

function [image_center_in, alpha, beta, gamma, nas, par, pal ] = ...
            calc_rot_angles( full_ref_local_file_name )

% Read in the local coordinates using the read function.
    [ref_locals] = read_ref_landmarks_file(full_ref_local_file_name) ;
    nas = ref_locals ( 1 , :) ;
    par = ref_locals ( 2 , :) ;
    pal = ref_locals ( 3 , :) ;

    % shift local reference points so par is at origin
    nas_or = nas - par ;
    pal_or = pal - par ;
    par_or = par - par ;
        
% Use PAR as reference; represented in rcp format for tformarray
% Is above statement rigt?
    image_center_in = par;
    
% Left Handed Analyze coordinate system
% Calculate each of the rotation angles
% recompute points after each rotation in order to get new angles
    
% rotate about z, positive rotation moves x into y    
    par_pal_xydist = sqrt((pal_or(1)-par_or(1))^2+(pal_or(2)-par_or(2))^2);
    par_nas_xydist = sqrt((nas_or(1)-par_or(1))^2+(nas_or(2)-par_or(2))^2);

    gamma = asind((par_or(2)-pal_or(2))/ par_pal_xydist) ; 
    
% Points (nas,pal) move after each rotation and angles must be
% recalculated.
    nas_gamma = asind(nas_or(2)/ par_nas_xydist);
    nas_rotz = nas_gamma + gamma ;
    
    pal_or(2) = sind(gamma-gamma)* par_pal_xydist;  
    pal_or(1) = cosd(gamma-gamma)* par_pal_xydist;
    nas_or(2) = sind(nas_rotz)* par_nas_xydist;  
    nas_or(1) = cosd(nas_rotz)* par_nas_xydist;
    
    % rotate about y, positive rotation moves z into x 
    par_pal_zxdist = sqrt((pal_or(1)-par_or(1))^2+(pal_or(3)-par_or(3))^2);
    par_nas_zxdist = sqrt((nas_or(1)-par_or(1))^2+(nas_or(3)-par_or(3))^2);
    
    beta = asind((-1*(par_or(3)-pal_or(3))/ par_pal_zxdist));
    
% Points (nas,pal) move after each rotation and angles must be
% recalculated. 
    nas_beta = asind(nas_or(3)/ par_nas_zxdist);
    nas_roty = nas_beta - beta ;
    
    pal_or(3) = sind(beta-beta)* par_pal_zxdist;  
    pal_or(1) = cosd(beta-beta)* par_pal_zxdist;
    nas_or(3) = sind(nas_roty)* par_nas_zxdist;  
    nas_or(1) = cosd(nas_roty)* par_nas_zxdist;
    
% rotate about x, positive rotation moves y into z    
    %par_pal_yzdist = sqrt((pal(2)-par(2))^2+(pal(3)-par(3))^2);
    par_nas_yzdist = sqrt((nas_or(2)-par_or(2))^2+(nas_or(3)-par_or(3))^2);
    
    alpha = asind((par_or(3)-nas_or(3))/ par_nas_yzdist);
     
% Points nas moves after rotation and angles must be
% recalculated.
    nas_or(3) = sind(alpha-alpha)* par_nas_yzdist;
    nas_or(2) = cosd(alpha-alpha)* par_nas_yzdist;
    
    nas = nas_or + par ;
    par = par_or + par ;
    pal = pal_or + par ;
       
end

