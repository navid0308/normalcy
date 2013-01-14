%**************************************************************************
%
%  PROGRAM TITLE      disp_coronal_movie.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       MAY 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%       disp_coronal_movie(g_3d, s_3d, b_3d_with_patch, pause_before, ...
%           pause_between) ;
%
%     where
%        g_3d            is an image volume containing the grayscale image
%                        from the CT scan. (double)
%        s_3d            is a binary-valued image volume containing the
%                        segmented head from the CT scan. s_3d is equal to
%                        0 outside the head and equal to 1 inside the head.
%                        (double)
%        b_3d_with_patch is a ternary-valued image volume equal to 0 for
%                        voxels not on the border of the segmented volume,
%                        equal to 1 for voxels on the border, and equal to
%                        2 for voxels on the border which are also within
%                        the patch of voxels to which the model sphere will
%                        be fitted. (double)
%        pause_before    is the number of seconds of delay between the time
%                        the blank figure is created and when the display
%                        of the movie begins. The purpose of this parameter
%                        is to allow the user to move the figure window and
%                        resize it during demonstrations of the software.
%                        Set this parameter to 0 for maximal speed.
%                        (double)
%        pause_between   is the number of seconds of delay between frames
%                        of the movie. Set this parameter to 0 for maximal
%                        speed. (double)
%
%  PROGRAM DESCRIPTION
%          disp_coronal_movie displays a three-panel movie of the coronal
%     slices derived from the CT images. The left-hand panel is the
%     gray-scale image, the middle panel is the segmentation of the head,
%     and the right-hand penel shows the border voxels in white and the
%     patch voxels in red.
%          A custom colormap is needed to properly display the red voxels
%     of the patch shown in the right-hand panel.
%          By default, Matlab displays row 1 of the image at the top. This
%     default has been overridden so that the image is displayed "correct
%     side up." The origin of the coordinate system is at the lower-left.
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

function [] = disp_coronal_movie(g_3d, s_3d, b_3d_with_patch, ...
    pause_before, pause_between)

    temp = size(g_3d) ;
    num_row = temp(1) ;
    
    % Create a custom colormap. This is necessary because we're using a
    % gray-scale colormap but want to display the patch voxels in red.

    gray_custom = gray(64) ;
    gray_custom(end,:) = [1 0 0] ;
    
    % Create a figure and pause to allow it to be moved and/or resized.

    figure(3) ;
    pause(pause_before) ;
    
    % Play a movie of successive coronal planes.
    
    for i = 1:num_row
        
        % Display the grayscale image on the left.

        subplot(1, 3, 1) ;
        imagesc(min(squeeze(g_3d(i ,:, :))', 1500)) ;
        lct = length(gray_custom) ;
        caxis_min = -1025 ;
        m = max(max(min(squeeze(g_3d(i, :, :))', 1500))) ;
        caxis_max = ceil((1 / (lct - 1)) * ((m * lct) - caxis_min)) ;
        caxis([caxis_min caxis_max]) ;
        colormap(gray_custom) ;
        h1 = gca ;
        set(h1, 'YDir', 'normal') ;
        axis equal ; axis tight ;
        xlabel('x (column)') ; ylabel('z (plane)') ;
        title(['Coronal plane ' num2str(i) ' of ' num2str(num_row)]) ;
        
        % Display the binary segmentation image in the middle.

        subplot(1, 3, 2) ;
        imagesc(squeeze(s_3d(i, :, :))') ;
        caxis([0 1.02]) ;
        colormap(gray_custom) ;
        h2 = gca ;
        set(h2, 'YDir', 'normal') ;
        axis equal ; axis tight ;
        xlabel('x (column)') ; ylabel('z (plane)') ;;
        title(['Coronal plane ' num2str(i) ' of ' num2str(num_row)]) ;
        
        % Display the border voxels in white and the patch voxels in red on
        % the right.

        subplot(1, 3, 3) ;
        imagesc(squeeze(b_3d_with_patch(i, :, :))') ;
        caxis([0 1.02]) ;
        colormap(gray_custom) ;
        h3 = gca ;
        set(h3, 'YDir', 'normal') ;
        axis equal ; axis tight ;
        xlabel('x (column)') ; ylabel('z (plane)') ;
        title(['Coronal plane ' num2str(i) ' of ' num2str(num_row)]) ;
        
        % Pause between frames of the movie (if pause_between is not equal
        % to 0).

        drawnow ;
        pause(pause_between) ;
        
    end

end