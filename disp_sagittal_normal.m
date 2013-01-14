%**************************************************************************
%
%  PROGRAM TITLE      disp_sagittal_normal.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  CALLING SYNTAX
%     Use the following syntax:
%       disp_sagittal_normal(g_3d, b_3d, p_3d, ref_point_rcp, ...
%           center_est_rcp, pause_before) ;
%
%     where
%        g_3d           is an image volume containing the grayscale image
%                       from the CT scan. (double)
%        b_3d           is a binary-valued image volume equal to 0 for
%                       voxels not on the border of the segmented volume,
%                       and equal to 1 for voxels on the border. (double)
%        p_3d           is a binary-valued image volume equal to 0 for
%                       voxels not in the patch of voxels to which the
%                       model sphere was fitted, and equal to 1 for voxels
%                       in the patch. (double)
%        ref_point_rcp  is a 1x3 row vector which is the reference point in
%                       rcp (row, column, plane) format. The reference
%                       point should be on the surface of the scalp and is
%                       the point at which the surface normal was
%                       calculated. (double)
%        center_est_rcp is a 1x3 row vector which is the center in rcp
%                       (row, column, plane) format of the best-fit sphere
%                       to the patch of voxels on the border. (double)
%        pause_before   is the number of seconds of delay between the time
%                       the blank figure is created and its display begins.
%                       The purpose of this parameter is to allow the user
%                       to move the figure window and resize it during
%                       demonstrations of the software. Set this parameter
%                       to 0 for maximal speed. (double)
%
%  PROGRAM DESCRIPTION
%          disp_sagittal_normal displays a three-panel figure of the
%     sagittal slice that contains the current reference point. The
%     left-hand panel is the gray-scale image, the middle panel shows the
%     border voxels of the segmentation of the left-hand panel, and the
%     right-hand panel shows the patch voxels. The patch voxels are a
%     subset of the border voxels which are used in the least-squares fit
%     of a sphere to find the surface normal at the reference point.
%          The reference point is shown as a red square and the estimated
%     center of the best-fit sphere is shown as a red circle in all three
%     panels. The circle and square are connected by a radial vector. This
%     radial vector should be perpendicular to the locus of patch voxels
%     shown in the right-hand panel.
%          By default, Matlab displays row 1 of the image at the top and
%     column 1 of the image at the left. Both of these defaults are
%     overridden so that the image is displayed "correct side up" and so
%     that the displayed head is facing the left, which is the standard for
%     sagittal images.
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

function [] = disp_sagittal_normal(g_3d, b_3d, p_3d, ref_point_rcp, ...
    center_est_rcp, pause_before)

    % Create a figure and pause to allow it to be moved and/or resized.
    
    figure(5) ;
    pause(pause_before) ;
    
    % Display the grayscale image on the left. Plot the reference point as
    % a red square, the projection into this plane of the center of the
    % best-fit sphere as a red circle, and draw a red line between them.

    subplot(1, 3, 1) ;
    hold off ;
    imagesc(squeeze(g_3d(:, ref_point_rcp(2), :))') ;
    colormap('gray') ;
    h1 = gca ;
    set(h1, 'XDir', 'reverse') ;
    set(h1, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('y (row)') ; ylabel('z (plane)') ;
    title(['Sagittal plane ' num2str(ref_point_rcp(2))]) ;
    hold on ;
    plot([center_est_rcp(1) ref_point_rcp(1)], ...
         [center_est_rcp(3) ref_point_rcp(3)], 'red') ;
    plot(ref_point_rcp(1), ref_point_rcp(3), 'rs') ;
    plot(center_est_rcp(1), center_est_rcp(3), 'ro') ;
    
    % Display the border voxels in white. Plot the reference point as a red
    % square, the projection into this plane of the center of the best-fit
    % sphere as a red circle, and draw a red line between them.

    subplot(1, 3, 2) ;
    hold off ;
    imagesc(squeeze(b_3d(:, ref_point_rcp(2), :))') ;
    colormap('gray') ;
    h2 = gca ;
    set(h2, 'XDir', 'reverse') ;
    set(h2, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('y (row)') ; ylabel('z (plane)') ;
    title(['Sagittal plane ' num2str(ref_point_rcp(2))]) ;
    hold on ;
    plot([center_est_rcp(1) ref_point_rcp(1)], ...
         [center_est_rcp(3) ref_point_rcp(3)], 'red') ;
    plot(ref_point_rcp(1), ref_point_rcp(3), 'rs') ;
    plot(center_est_rcp(1), center_est_rcp(3), 'ro') ;

    % Display the patch voxels in white. Plot the reference point as a red
    % square, the projection into this plane of the center of the best-fit
    % sphere as a red circle, and draw a red line between them.
    
    subplot(1, 3, 3) ;
    hold off ;
    imagesc(squeeze(p_3d(:, ref_point_rcp(2), :))') ;
    colormap('gray') ;
    h3 = gca ;
    set(h3, 'XDir', 'reverse') ;
    set(h3, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('y (row)') ; ylabel('z (plane)') ;
    title(['Sagittal plane ' num2str(ref_point_rcp(2))]) ;
    hold on ;
    plot([center_est_rcp(1) ref_point_rcp(1)], ...
         [center_est_rcp(3) ref_point_rcp(3)], 'red') ;
    plot(ref_point_rcp(1), ref_point_rcp(3), 'rs') ;
    plot(center_est_rcp(1), center_est_rcp(3), 'ro') ;

end