%**************************************************************************
%
%  PROGRAM TITLE     calc_circumf_angles.m
%
%  WRITTEN BY        Gregory G. Reiker and Kirk Smith
%  DATE WRITTEN      September 14, 2011 
%  WRITTEN FOR       Pediatric head modeling project
%
%  REVISIONS BY      
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [image_center_in, alpha, beta, gamma, fp, op ] = ...
%           calc_circumf_angles(fp, op ) ;
%
%     where
%       image_center_in row, column, plane of image center.
%       alpha           angle of rotation about x.
%       beta            angle of rotation about y.
%       gamma           angle of rotation about z.
%       fp              output frontal pole.
%       op              output occipital pole.
%       fp              input frontal pole.
%       op              input occipital pole.
%
%  PROGRAM DESCRIPTION
%      This program will be called from within Calc_braincase since the
%   segmented skull and bounding box coordinates are already calculated
%   there. The maximum y coordinate (Analyze coordinate system so may need
%   to adjust based on r,c,p of Matlab) is the frontal pole. The minimum
%   Analyze y coordinate is the occipital pole. This should not be called
%   until after the volume has been aligned to the nas, par, pal plane
%   which is done in calc_rot_angles.
%       Calculate a Euler angles to align the frontal pole and occipital poles
%   automatically. When applying the rotations, order of execution matters.
%   Make sure to follow it. First step is to calculate the rotations in
%   degrees. We will use the occipital pole (op) as the reference,
%   therefore, par_or will be op and nas_or will be fp (frontal pole). The
%   only rotation performed will be about x, so set beta and gamma to 0. 
%       Based on calc_rot_angles.    
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
%       1.0     September 14, 2011   Initial release.
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


function [image_center_in, alpha, beta, gamma, fp, op ] = ...
            calc_circumf_angles(fp, op )

% Read in the local coordinates using the read function.
    nas = fp ;
    par = op ;

    % shift local reference points so par is at origin
    nas_or = nas - par ;
    par_or = par - par ;
        
% Use PAR as reference; represented in rcp format for tformarray
    image_center_in = par;
    
% Left Handed Analyze coordinate system
% Calculate each of the rotation angles
% recompute points after each rotation in order to get new angles
    
    gamma = 0 ; 
    
    beta = 0;
    
% rotate about x, positive rotation moves y into z    
    par_nas_yzdist = sqrt((nas_or(2)-par_or(2))^2+(nas_or(3)-par_or(3))^2);
    
    alpha = asind((par_or(3)-nas_or(3))/ par_nas_yzdist);
     
% Points nas moves after rotation and angles must be
% recalculated.
    nas_or(3) = sind(alpha-alpha)* par_nas_yzdist;
    nas_or(2) = cosd(alpha-alpha)* par_nas_yzdist;
    
    fp = nas_or + par ;
    op = par_or + par ;
    
       
end

