%**************************************************************************
%
%  PROGRAM TITLE      fit_function_sphere.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%       [sse] = fit_function_sphere(parms_in, row_data, col_data, pln_data) ;
%
%     where
%        sse      is the "sum of squared errors," which is the objective
%                 function used when fitting a sphere to the patch of
%                 voxels on the border of the segmentation of the head. The
%                 sum of squared errors is defined in spherical
%                 coordinates, and is equal to the sum of the squares of
%                 the differences between the radius of the data points (on
%                 the patch) and the radius of the circle that is being
%                 fit. (double)
%        parms_in is a 4x1 column vector of parameters that defines the
%                 sphere:
%                    parms_in(1) is the (possibly fractional) row
%                                coordinate of the center of the circle;
%                    parms_in(2) is the (possibly fractional) column
%                                coordinate of the center of the circle;
%                    parms_in(3) is the (possibly fractional) plane
%                                coordinate of the center of the circle;
%                    parms_in(4) is the radius of the circle in (row,
%                                column, plane) units of length.
%        row_data is a vector with n elements, where n is the number of
%                 voxels in the patch. Each element of row_data is the row
%                 index of a voxel in the patch.
%        col_data is a vector with n elements, where n is the number of
%                 voxels in the patch. Each element of col_data is the
%                 column index of a voxel in the patch.
%        pln_data is a vector with n elements, where n is the number of
%                 voxels in the patch. Each element of pln_data is the
%                 plane index of a voxel in the patch.
%
%  PROGRAM DESCRIPTION
%          fit_function_sphere is the objective function that is minimized
%     to find the sphere that best fits the patch of voxels that have been
%     defined on the border of the segmentation of the head. The objective
%     function is a "sum of squared errors" measure implemented in
%     spherical coordinates.
%          This function is called repeatedly by the Optimization Toolbox
%     function called fminsearch.
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
%  REVISION HISTORY
%     Version      Date                          Comment
%     -------   -----------   ---------------------------------------------
%       1.0     12 May 2011   Initial release.
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

function [sse] = ...
    fit_function_sphere(parms_in, row_data, col_data, pln_data)

    row_center = parms_in(1) ;
    col_center = parms_in(2) ;
    pln_center = parms_in(3) ;
    radius     = parms_in(4) ;

    sse = sum((radius - sqrt(((row_data - row_center).^2) ...
                           + ((col_data - col_center).^2) ...
                           + ((pln_data - pln_center).^2))).^2) ;

end