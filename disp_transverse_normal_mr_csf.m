%**************************************************************************
%
%  PROGRAM TITLE      disp_transverse_normal_mr_csf.m
%
%  WRITTEN BY         David G. Politte
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY       Gregory G. Reiker
%  DATE MODIFIED      July 28, 2011
%  DATE MODIFIED      November 10, 2011
%
%  CALLING SYNTAX
%     Use the following syntax:
%       disp_transverse_normal_mr(mr_3d, g_3d, b_3d, p_3d, ref_point_rcp, ...
%           center_est_rcp, pause_before) ;
%
%     where
%        mr_3d          is an image volume containing the grayscale image
%                       from the MR scan. (double)
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
%        bi_edge_rcp    is the inner bone edge rcp.
%        csf_edge_rcp   is the CSF edge rcp.
%        pause_before   is the number of seconds of delay between the time
%                       the blank figure is created and its display begins.
%                       The purpose of this parameter is to allow the user
%                       to move the figure window and resize it during
%                       demonstrations of the software. Set this parameter
%                       to 0 for maximal speed. (double)
%
%  PROGRAM DESCRIPTION
%          disp_transverse_normal displays a 4-panel figure of the
%     transverse slice that contains the current reference point. The
%     upper left-hand panel is the MR gray-scale image, the right is the CR
%     gray-scale image,the lower left panel shows the
%     border voxels of the segmentation of the CR panel, and the
%     lower right panel shows the patch voxels. The patch voxels are a
%     subset of the border voxels which are used in the least-squares fit
%     of a sphere to find the surface normal at the reference point.
%          The reference point is shown as a red square and the estimated
%     center of the best-fit sphere is shown as a red circle in all three
%     panels. The circle and square are connected by a radial vector. This
%     radial vector should be perpendicular to the locus of patch voxels
%     shown in the right-hand panel.
%          By default, Matlab displays row 1 of the image at the top. This
%     default has been overridden so that row 1 is displayed at the bottom.
%     The origin of the coordinate system is at the lower-left.
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
%       1.1     27 July 2011  Added MR plot.
%       1.2     November 10, 2011  Added points for bone and CSF edges.
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

function [] = disp_transverse_normal_mr_csf(mr_3d, g_3d, b_3d, p_3d, ref_point_rcp, ...
    center_est_rcp, bi_edge_rcp, csf_edge_rcp, pause_before)

    % Create a figure and pause to allow it to be moved and/or resized.

    figure(4) ;
    pause(pause_before) ;
    
   % Display the MR grayscale image on the left. Plot the reference point as
    % a red square, the projection into this plane of the center of the
    % best-fit sphere as a red circle, and draw a red line between them.

    subplot(2, 2, 1) ;
    hold off ;
    imagesc(mr_3d(: ,:, ref_point_rcp(3))) ;
    colormap('gray') ;
    h1 = gca ;
    set(h1, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('x (column)') ; ylabel('y (row)') ;
    title(['Transverse plane ' num2str(ref_point_rcp(3))]) ;
    hold on ;
    plot([center_est_rcp(2) ref_point_rcp(2)], ...
         [center_est_rcp(1) ref_point_rcp(1)], 'red') ;
    plot(ref_point_rcp(2), ref_point_rcp(1), 'rs') ;
    plot(center_est_rcp(2), center_est_rcp(1), 'ro') ;
    plot(bi_edge_rcp(2), bi_edge_rcp(1), 'gx');
    plot(csf_edge_rcp(2), csf_edge_rcp(1), 'yx');
    
    % Display the CT grayscale image next. Plot the reference point as
    % a red square, the projection into this plane of the center of the
    % best-fit sphere as a red circle, and draw a red line between them.

    subplot(2, 2, 2) ;
    hold off ;
    imagesc(g_3d(: ,:, ref_point_rcp(3))) ;
    colormap('gray') ;
    h1 = gca ;
    set(h1, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('x (column)') ; ylabel('y (row)') ;
    title(['Transverse plane ' num2str(ref_point_rcp(3))]) ;
    hold on ;
    plot([center_est_rcp(2) ref_point_rcp(2)], ...
         [center_est_rcp(1) ref_point_rcp(1)], 'red') ;
    plot(ref_point_rcp(2), ref_point_rcp(1), 'rs') ;
    plot(center_est_rcp(2), center_est_rcp(1), 'ro') ;
    plot(bi_edge_rcp(2), bi_edge_rcp(1), 'gx');
    plot(csf_edge_rcp(2), csf_edge_rcp(1), 'yx');
    
    % Display the border voxels in white. Plot the reference point as a red
    % square, the projection into this plane of the center of the best-fit
    % sphere as a red circle, and draw a red line between them.

    subplot(2, 2, 3) ;
    hold off ;
    imagesc(b_3d(:, :, ref_point_rcp(3))) ;
    colormap('gray') ;
    h2 = gca;
    set(h2, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('x (column)') ; ylabel('y (row)') ;
    title(['Transverse plane ' num2str(ref_point_rcp(3))]) ;
    hold on ;
    plot([center_est_rcp(2) ref_point_rcp(2)], ...
         [center_est_rcp(1) ref_point_rcp(1)], 'red') ;
    plot(ref_point_rcp(2), ref_point_rcp(1), 'rs') ;
    plot(center_est_rcp(2), center_est_rcp(1), 'ro') ;
    
    % Display the patch voxels in white. Plot the reference point as a red
    % square, the projection into this plane of the center of the best-fit
    % sphere as a red circle, and draw a red line between them.

    subplot(2, 2, 4) ;
    hold off ;
    imagesc(p_3d(:, :, ref_point_rcp(3))) ;
    colormap('gray') ;
    h3 = gca ;
    set(h3, 'YDir', 'normal') ;
    axis equal ; axis tight ;
    xlabel('x (column)') ; ylabel('y (row)') ;
    title(['Transverse plane ' num2str(ref_point_rcp(3))]) ;
    hold on ;
    plot([center_est_rcp(2) ref_point_rcp(2)], ...
         [center_est_rcp(1) ref_point_rcp(1)], 'red') ;
    plot(ref_point_rcp(2), ref_point_rcp(1), 'rs') ;
    plot(center_est_rcp(2), center_est_rcp(1), 'ro') ;
  
end