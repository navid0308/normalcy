%**************************************************************************
%
%  PROGRAM TITLE      calc_landmarks.m
%
%  WRITTEN BY         Kirk E. Smith
%  DATE WRITTEN       September 20, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY      
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [ ref_point_mat_xyz, num_ref_points, local_origin ] = ...
%           calc_landmarks( ref_landmark_file_name, list_border_bc_rcp, ...
%                    min_rcp_of_data, max_rcp_of_data, nas, par, pal) ;
%
%     where
%       ref_point_mat_xyz   matrix of x,y,z of reference points
%       num_ref_points      number of reference points
%       local_origin        x,y,z of local origin
%       ref_landmark_file_name name of the text file that contains
%                             the (azimuth_deg,elevation_deg,radius_mm)
%                             values of the reference points that are used
%                             to define the landmarks. default r=1.
%       list_border_bc_rcp  list of border before cleanup r,c,p
%       min_rcp_of_data     minimum row, column, plane of data
%       max_rcp_of_data     maximum row, column, plane of data
%       nas                 NAS
%       par                 PAR
%       pal                 PAL
%
%  PROGRAM DESCRIPTION
%    Function calculates the landmarks. It reads the nas, par, and pal 
% coordinates from the output of calc_rot_angles to reflect the 
% auto-alignment process of nas, par, and pal. Added nas, par, pal as inputs.
%   Min distance point to line - Dave Politte.
%   Beware as rcp means row, column, plane and not x,y,z.  
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
%       read_ref_landmarks_file
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     September 20, 2011   Initial release.
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


function [ ref_point_mat_xyz, num_ref_points, local_origin ] = ...
    calc_landmarks( ref_landmark_file_name, list_border_bc_rcp, ...
                    min_rcp_of_data, max_rcp_of_data, nas, par, pal)

 % read_ref_landmarks_file is a function used to read in x,y,z coordinates
% from a text file. ref_local_file_name is the local coordinate system
% defined by NAS, PAR, PAL. ref_landmark_file_name contains the azimuth,
% and elevation angles in degrees that define the 10-20 array landmarks.
% The 3rd column is set as 1 by default so that the function
% read_ref_landmarks_file can be used without modification, but the 3rd
% column (r) is calculated from the local coordinates.    
%
% These coordinates are no longer appropriate. The nas, par, and pal
% landmarks now need to come from calc_roataion_angles. Simply pass the
% values in as they are returned to normalcy from calc_rotation_angles.
% [ref_locals] = read_ref_landmarks_file(ref_local_file_name) ;
% nas = ref_locals ( 1 , :) ;
% par = ref_locals ( 2 , :) ;
% pal = ref_locals ( 3 , :) ;
[ref_landmarks] = read_ref_landmarks_file(ref_landmark_file_name) ;

azimuth_deg = ref_landmarks( : , 1) ;
elevation_deg = ref_landmarks( : , 2) ;

% convert to radians
azimuth_rad = azimuth_deg * pi / 180 ;
elevation_rad = elevation_deg * pi / 180 ;

% calculate local origin and radius from ref_locals
local_r = (pal(1) - par(1)) / 2 ;
local_x = (par(1) + pal(1)) / 2 ;
%local_y = (par(2) + pal(2)) / 2 ;
local_y = (min_rcp_of_data(1) + max_rcp_of_data(1)) / 2 ;
local_z = par(3) ;
local_origin = [local_x local_y local_z] ;

% convert azimuth, elevation, radius to cartesian coordinates
cart_ref_landmarks = ...
    [local_x + (local_r * cos(elevation_rad) .* cos(azimuth_rad))...
    local_y + (local_r * cos(elevation_rad) .* sin(azimuth_rad))...
    local_z + (local_r * sin(elevation_rad))] ;
 
% calculate minimum distance between border point and line.
% For each border point calculate min dist to line. Line is established
% using local_origin and cart_ref_landmarks.
%
%
% QUESTION: What is the shortest distance between a point (x,y,z)
% and the line through (x1,y1,z1) and (x2,y2,z2)?
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ANSWER: Let d0 be the desired shortest distance. Then d0 is
% given by
% 
%      d0 = sqrt[(s(x2 - x1) + x1 - x)^2
%               +(s(y2 - y1) + y1 - y)^2
%               +(s(z2 - z1) + z1 - z)^2],
% with
% 
%     (x2 - x1)(x - x1) + (y2 - y1)(y - y1) + (z2 - z1)(z - z1)
% s = --------------------------------------------------------- .
%            (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2
%
% x = border_x (border pt)
% y = border_y (border pt)
% z = border_z (border pt)
% x1 = local_x
% y1 = local_y
% z1 = local_z
% x2 = cart_ref_landmarks (x value)
% y2 = cart_ref_landmarks (y value)
% z2 = cart_ref_landmarks (z value)
%
border_x = list_border_bc_rcp( : ,2) ;
border_y = list_border_bc_rcp( : ,1) ;
border_z = list_border_bc_rcp( : ,3) ;
landmark_x = cart_ref_landmarks( : ,1) ;
landmark_y = cart_ref_landmarks( : ,2) ;
landmark_z = cart_ref_landmarks( : ,3) ;
x1 = local_x ;
y1 = local_y ;
z1 = local_z ;
x2 = landmark_x ;
y2 = landmark_y ;
z2 = landmark_z ;
% x = border_x ;
% y = border_y ;
% z = border_z ;

% Find the border voxels excluding the ones less than local_z
% Added +4 (+50 temp) to allow for slight offsets of local z and first slice
[index] = find(border_z > local_z+4) ;
border_x = border_x(index) ;
border_y = border_y(index) ;
border_z = border_z(index) ;
list_border_bc_rcp = [border_x border_y border_z] ;

x = border_x ;
y = border_y ;
z = border_z ;

temp0 = size(list_border_bc_rcp) ;
row_tot = temp0(1) ;
temp1 = size(cart_ref_landmarks) ;
land_tot = temp1(1) ;
pt2line_dist = zeros(row_tot, 1) ;
ref_point_mat_xyz = zeros(land_tot, 3) ;
for land = 1:land_tot
    for row = 1:row_tot
     %do the main loop here
     %   border_pt_xyz = [list_border_bc_rcp(2) list_border_bc_rcp(1)...
     %       list_border_bc_rcp(3)] ;
      point_line_offset = ((x2(land) - x1) .* (x(row) - x1)...
          + (y2(land) - y1) .* (y(row) - y1)...
          + (z2(land) - z1) .* (z(row) - z1))...
          / ((x2(land) - x1).^2 + (y2(land) - y1).^2 ...
          + (z2(land) - z1).^2) ;
    s = point_line_offset ;
    pt2line_dist(row) = sqrt((s .* (x2(land) - x1) + x1 - x(row)).^2 ...
        +(s .* (y2(land) - y1) + y1 - y(row)).^2 ...
        +(s .* (z2(land) - z1) + z1 - z(row)).^2) ;

    end
    [short_dist,short_index] = min(pt2line_dist) ;
    ref_point_mat_xyz(land, : ) = [border_x(short_index)...
        border_y(short_index)...
        border_z(short_index)] ;
end
temp = size(ref_point_mat_xyz) ;
    num_ref_points = temp(1) ;
    if (num_ref_points <= 0)
        disp('ERROR:  NO REFERENCE POINTS SPECIFIED. EXITING!') ;
        return ;
    end