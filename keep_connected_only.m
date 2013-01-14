%**************************************************************************
%
%  PROGRAM TITLE      keep_connected_only.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [img_out_3d, rcp_list_out] = ...
%           keep_connected_only(img_in_3d, rcp_seed) ;
%
%     where
%        img_out_3d   is a binary valued image volume equal to 0 in some
%                     voxels and equal to 1 in voxels which were equal to 1
%                     in img_in_3d and which are connected to the voxel
%                     specified in rcp_seed. (double)
%        rcp_list_out is an nx3 matrix that contains the (row, column,
%                     plane) coordcinates of all of the voxels in
%                     img_out_3d that are equal to 1. Here, n is the number
%                     of voxels equal to 1. The two outputs of this Matlab
%                     function, img_out_3d and rcp_list_out, contain
%                     exactly the same information but in different
%                     formats. (double)
%        img_in_3d    is a binary valued image volume equal to 0 in some
%                     voxels and equal to 1 in others. (double)
%        rcp_seed     is a 1x3 row vector containing the row, column, and
%                     plane indices of a voxel that is within the set of
%                     voxels in img_in_3d that have value equal to 1.
%                     (double)
%
%  PROGRAM DESCRIPTION
%          keep_connected_only is used to "clean up" the border of the
%     segmentation or the patch of voxels to which a sphere will be fit.
%          Given a binary valued input image volume, img_in_3d, and the
%     indices of a voxel within img_in_3d that has value equal to 1,
%     rcp_seed, this program outputs a binary valued image volume,
%     img_out_3d, that is equal to 1 in all the voxels that are both (1)
%     equal to 1 in img_in_3d, and (2) are connected to the voxel specified
%     in rcp_seed. By "connected," we mean that the voxel has either
%     north-south, east-west, or front-back adjacency to another voxel that
%     is equal to 1 in img_in_3d, and that this type of chained adjacency
%     leads in one or more steps to the seed voxel.
%          The segmentations typically contain some "noisy" voxels such
%     that the border finding algorithm includes them in the border
%     although they are not truly on the exterior of the head. This region-
%     growing algorithm is applied twice, once when defining the overall
%     border, and once when defining the patch, to clean up and reduce the
%     set of voxels on the border or the patch, respectively.
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

function [img_out_3d, rcp_list_out] = keep_connected_only(img_in_3d, rcp_seed)

    % Create the output image.
    
    img_out_3d = zeros(size(img_in_3d)) ;
  
    % Initialize arrays and variables needed by the region growing
    % algorithm. Start with the seed voxel.
    
    new_grow_list = rcp_seed ;
    img_out_3d(new_grow_list(1), new_grow_list(2), new_grow_list(3)) = 1 ;

    temp = size(new_grow_list) ;
    num_grow = temp(1) ;

    % Keep adding voxels to the output list until the size of the list
    % doesn't change any more, which means we're finished. The vector
    % grow_list contains the set of all voxels that were added to the
    % output the last time through the loop.
    
    while (num_grow ~= 0)
        
        grow_list = new_grow_list ;
        new_grow_list = [] ;
        
        for m = 1:num_grow
            
           i = grow_list(m, 1) ;
           j = grow_list(m, 2) ;
           k = grow_list(m, 3) ;
           
           % Add voxels to the output list if and only if they are adjacent
           % to a voxel that is already in the output list and is in the
           % current "grow list." We need to check all 6 directions for
           % adjacency.
           
           if (and((img_in_3d(i - 1, j, k) == 1), ...
                   (img_out_3d(i - 1, j, k) ~= 1)))
               img_out_3d(i - 1, j, k) = 1 ;
               new_grow_list = [new_grow_list ; i-1 j k] ;
           end
           if (and((img_in_3d(i + 1, j, k) == 1), ...
                   (img_out_3d(i + 1, j, k) ~= 1)))
               img_out_3d(i + 1, j, k) = 1 ;
               new_grow_list = [new_grow_list ; i+1 j k] ;
           end
           if (and((img_in_3d(i, j - 1, k) == 1), ...
                   (img_out_3d(i, j - 1, k) ~= 1)))
               img_out_3d(i, j - 1, k) = 1 ;
               new_grow_list = [new_grow_list ; i j-1 k] ;
           end
           if (and((img_in_3d(i, j + 1, k) == 1), ...
                   (img_out_3d(i, j + 1, k) ~= 1)))
               img_out_3d(i, j + 1, k) = 1 ;
               new_grow_list = [new_grow_list ; i j+1 k] ;
           end
           if (and((img_in_3d(i, j, k - 1) == 1), ...
                   (img_out_3d(i, j, k - 1) ~= 1)))
               img_out_3d(i, j, k - 1) = 1 ;
               new_grow_list = [new_grow_list ; i j k-1] ;
           end
           if (and((img_in_3d(i, j, k + 1) == 1), ...
                   (img_out_3d(i, j, k + 1) ~= 1)))
               img_out_3d(i, j, k + 1) = 1 ;
               new_grow_list = [new_grow_list ; i j k+1] ;
           end
           
        end
        
        temp = size(new_grow_list) ;
        num_grow = temp(1) ;
        
    end
    
    % Build lists of the rows, columns, and planes of voxels that are equal
    % to 1 in the input image and are connected to the seed voxel.
    
    [i, j, k] = ind2sub(size(img_out_3d), find(img_out_3d > 0)) ;
    rcp_list_out = [i j k] ;

end