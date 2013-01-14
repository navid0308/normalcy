%**************************************************************************
%
%  PROGRAM TITLE      calc_border_bbox.m
%
%  WRITTEN BY         Kirk E. Smith
%  DATE WRITTEN       November 3, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%       [ min_rcp_of_data, max_rcp_of_data, list_border_bc_rcp, b_3d_bc  ] = ...
%          calc_border_bbox( s_3d );
%
%     where
%       min_rcp_of_data    minimum row, column, plane of data
%       max_rcp_of_data    maximum row, column, plane of data
%       list_border_bc_rcp list of border before cleanup r,c,p
%       b_3d_bc             border box 3D image before cleanup
%       s_3d                segmented image volume
%
%  PROGRAM DESCRIPTION
%   This function calculates the border bounding box based on code from normalcy.m .          
%
%  FILES
%     standard input - not used
%
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
%**************************************************************************

function [ min_rcp_of_data, max_rcp_of_data, list_border_bc_rcp, b_3d_bc  ] = ...
    calc_border_bbox( s_3d )

        disp('Calculating a bounding box for the segmentation ...') ;
        [temp1, temp2, temp3] = ind2sub(size(s_3d), find(s_3d > 0)) ;
        list_bin_rcp = [temp1 temp2 temp3] ;
        min_rcp_of_data = [min(list_bin_rcp(:,1)) ...
                       min(list_bin_rcp(:,2)) ...
                       min(list_bin_rcp(:,3))] ;
        max_rcp_of_data = [max(list_bin_rcp(:,1)) ...
                       max(list_bin_rcp(:,2)) ...
                       max(list_bin_rcp(:,3))] ;

        % Find all the border voxels.
        disp('Finding border voxels ...') ;
        kern = ones(3,3,3) ;
        % Below, "bc" means "before cleanup".
        b_3d_bc = sign((1 - s_3d) .* convn(s_3d, kern, 'same')) ;
        [temp1, temp2, temp3] = ind2sub(size(b_3d_bc), find(b_3d_bc > 0)) ;
        list_border_bc_rcp = [temp1 temp2 temp3] ;
 
        clear temp1 temp2 temp3 ;
end

