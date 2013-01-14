%**************************************************************************
%
%  PROGRAM TITLE      normalcy_master.m
%
%  WRITTEN BY         David G. Politte and Kirk E. Smith
%  DATE WRITTEN       May 12, 2011 
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY       Kirk Smith
%                     Gregory G. Reiker
%  DATE MODIFIED      2011, 2012
%
%  CALLING SYNTAX
%     Use the following syntax:
%        normalcy_master ;
%
%  PROGRAM DESCRIPTION
%       This script function specifies all of the variables needed to call
%   normalcy.m for a particular CT/MR case.  
%       It needs to be modified for each case run. 
%
%  FILES
%     standard input - not used
%     standard output - not used
%
%  DEPENDENCIES
%         MATLAB     (win64)                Version 7.12.0.635 (R2011a)
%         Image Acquisition Toolbox         Version 4.1        (R2011a)
%         Image Processing Toolbox          Version 7.2        (R2011a)
%         Optimization Toolbox              Version 6.0        (R2011a)
%         Parallel Computing Toolbox        Version 5.1        (R2011a)
%         Signal Processing Toolbox         Version 6.15       (R2011a)
%         Statistics Toolbox                Version 7.5    
%
%     Code dependencies are indicated by the level of indent:
%     normalcy_master
%       normalcy
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     May 12, 2011       Initial release.
%       1.1     November 1, 2011   Input variables for normalcy call changed.
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



inputs_directory = '' ;
outputs_directory = '' ;
% case1
% this was for manual entry of points: ref_points_file_name = 'experiment_y.log' ;
gray_file_name_root = '' ;
ref_local_file_name = '' ;
ref_landmark_file_name = '' ;
out_file_name_root = '' ;
mr_file_name_root = ''; 
mr_ref_landmark_file_name = 't' ; % 
r_patch_mm = 20 ;
display_verbose = 0 ; % If 0, don't display 3 movies, if 1, do display 3 movies
pause_before = 5 ;
pause_between = 0 ;
normalcy(inputs_directory, outputs_directory, ...
    gray_file_name_root, ref_local_file_name, ...
    ref_landmark_file_name, out_file_name_root, ...
    mr_file_name_root, ...
    mr_ref_landmark_file_name, r_patch_mm, ...
    display_verbose, pause_before, pause_between) ;
