%**************************************************************************
%
%  PROGRAM TITLE      read_ref_landmarks_file.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [ref_point_mat_xyz] = read_ref_points_file(file_name) ;
%
%     where
%        ref_point_mat_xyz is an nx3 matrix of the (x,y,z) coordinates of
%                          the reference points on the scalp where we wish
%                          to compute the direction vectors of the surface
%                          normals. Here, n is the number of reference
%                          points. (double)
%        file_name         is the name of the input file that contains the
%                          (x,y,z) coordinates of the reference points.
%                          This function is intended to read a log file
%                          created by Analyze software; the first 7 lines
%                          of text are header information that is skipped.
%                          (char)
%
%  PROGRAM DESCRIPTION
%          read_ref_points_file reads the (x,y,z) reference points from an
%     Analyze log file. The reference points are the points where we wish
%     to compute surface normals.
%          This function echos the file to the Matlab console.
%
%  FILES
%     standard input - not used
%     standard output - The input file, file_name, is echoed to the Matlab
%        console window.
%
%  DEPENDENCIES
%         MATLAB     (win64)                Version 7.12.0.635 (R2011a)
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

function [ref_point_mat_xyz] = read_ref_landmarks_file(file_name)

    % Open the input file with read-only privileges.
    
    fid = fopen(file_name, 'r') ;

    % Initialize variables, including the output array.
    
    num_points = 0 ;
    tline = 0 ;
    ref_point_mat_xyz = [] ;
    temp_mat = zeros(1, 3) ;

%     Skip over the first 7 header lines; they don't contain reference points.

    for i = 1:7
        tline = fgetl(fid) ;
        if (tline == -1)
            break ;
      end
    end

    % Read reference points until an end-of-file is encountered.

    while (tline ~= -1)
        tline = fgetl(fid) ;

        if (tline ~= -1)
            num_points = num_points + 1 ;
            temp_mat(num_points, :) = sscanf(tline, '%f', 3) ;
        end
    end
    
    if (num_points > 0)
       ref_point_mat_xyz = temp_mat ;
    end

    fclose(fid) ;
      
end