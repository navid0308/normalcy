%**************************************************************************
%
%  PROGRAM TITLE      normalcy.m
%
%  WRITTEN BY         David G. Politte and Kirk E. Smtih
%  DATE WRITTEN       May 12, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS          May 18, 2011; Kirk Smith 
%  REVISIONS          July 7, 2011; Gregory G. Reiker
%  REVISIONS          July 21, 2011; Gregory G. Reiker
%  REVISIONS          August 10, 2011; Kirk Smith 
%  REVISIONS          August 4-25, 2011; Gregory G. Reiker
%  REVISIONS          September 1-29, 2011; Gregory G. Reiker
%  REVISIONS          October 1-30, 2011; Kirk Smith
%  REVISIONS          October 6-27, 2011; Gregory G. Reiker
%  REVISIONS          November 1-30, 2011; Kirk Smith
%  REVISIONS          November 1-30, 2011; Gregory G. Reiker
%  REVISIONS          December 1-30, 2011; Kirk Smith
%  REVISIONS          December 1-30, 2011; Gregory G. Reiker
%
%  CALLING SYNTAX
%     Use the following syntax:
%      normalcy(inputs_directory, outputs_directory, ...
%        gray_file_name_root, ref_local_file_name, ...
%        ref_landmark_file_name, out_file_name_root, ...
%        mr_file_name_root, ...
%        mr_ref_landmark_file_name, r_patch_mm, ...
%        display_verbose, pause_before, pause_between) ;
%
%     where
%        inputs_directory     is the name of the directory (folder) where
%                             the inputs are located. The inputs are: (1)
%                             the header file of the CT grayscale image
%                             volume, (2) the image file of the CT grayscale
%                             image volume, (3) the header file of the
%                             MR grayscale image volume, (4)
%                             the image file of the MR grayscale image
%                             volume.
%                             The image volumes (items 1-4) are in Analyze 7.5
%                             format. (char)
%                             Revision (5-18-2011): 
%                             6) a text file that contains the (x,y,z)
%                             coordinates of the reference points that 
%                             define the local coordinate system (NAS,
%                             PAR, PAL)  
%                             7) a text file that contains azimuth and
%                             elevation angles to define landmarks where
%                             surface normals are calculated.                               
%        outputs_directory    is the name of the directory (folder) where
%                             the outputs will be located. The outputs are:
%                             (1) a ".mat" file that contains the contents
%                             of important variables, vectors, and
%                             matrices, (2) a ".csv" (comma separated
%                             values) list that contains the contents of
%                             important variables, vectors, and matrices
%                             and can be opened in Microsoft Excel. Items
%                             (1) and (2) are cumulative outputs for all of
%                             the reference points that were specified. In
%                             addition, for each reference point, the
%                             following outputs are created: (3) a ".png"
%                             graphic file of Figure 4, which is a
%                             three-panel representation of the transverse
%                             image slice of the reference point, (4) a
%                             ".png" graphic file of Figure 5, which is a
%                             three-panel representation of the sagittal
%                             image slice of the reference point, (5) a
%                             ".png" graphic file of Figure 6, which is a
%                             three-panel representation of the coronal
%                             image slice of the reference point, and (6) a
%                             ".png" file of the graph of the interpolated
%                             image values of the grayscale CT image along
%                             the computed surface normal. (char)
%        gray_file_name_root  is the root name of the grayscale CT image.
%                             The names of the header and image files are
%                             given by gray_file_name_root with ".hdr" or
%                             ".img" appended, respectively. (char) 
%        ref_local_file_name is the name of the text file that
%                             contains the (x,y,z) coordinates of the
%                             CT reference points (NAS,PAR,PAL) that define a
%                             local coordinate system.                             
%        ref_landmark_file_name is the name of the text file that contains
%                             the (azimuth_deg,elevation_deg,radius_mm)
%                             values of the reference points that are used
%                             to define the landmarks. default r=1.                            
%        out_file_name_root   is the root name of the output files that are
%                             generated. See "outputs_directory" above for
%                             a list of output files. (char)
%        mr_file_name_root    is the root name of the MR image.
%                             The names of the header and image files are
%                             given by mr_file_name_root with ".hdr" or
%                             ".img" appended, respectively. (char) 
%        mr_ref_landmark_file_name is the name of the text file that contains
%                             the (azimuth_deg,elevation_deg,radius_mm)
%                             values of the reference points that are used
%                             to define the landmarks for MR/CT registration. 
%                             default r=1.                            
%        r_patch_mm           is the radius in mm of the patch of voxels
%                             that will be included in the fit to a sphere.
%                             All voxels on the surface border within
%                             r_patch_mm of the reference point are
%                             included in the patch. (double)
%        display_verbose      indicates the degree of verbosity in the
%                             displayed results. If display_verbose equals
%                             1, three movies of the transverse, sagittal,
%                             and coronal planes are displayed. If
%                             display_verbose equals 0, these movies are
%                             not shown. (double)
%        pause_before         is the number of seconds between the time a
%                             figure window appears and the contents of the
%                             figure appear. This gives the user time to
%                             move and/or resize the window while giving
%                             demonstrations of the program. For maximal
%                             speed, set pause_before to 0 seconds.
%                             (double)
%        pause_between        is the number of seconds between frames of
%                             the transverse, sagittal, and coronal movies.
%                             For maximal speed, set pause_between to 0
%                             seconds. This parameter has no effect if
%                             display_verbose is 0. (double)
%
%  PROGRAM DESCRIPTION
%          normalcy is the main Matlab function for fitting spheres to
%     convex patches on the surface of the scalp. It also computes the
%     values of a grayscale CT image along the surface normals.
%          Also, it calculates braincase segmentation, accecpts MR data 
%     on which analyses are performed, performs 3D rotation and segmentation of
%     original grayscale CT image, calculates CSF measurement from MR data,
%     and performs 3D registration of MR image with CT image volume.
%
%  FILES
%     Inputs (located in "inputs_directory"):
%        grayscale CT image header (Analyze 7.5 format)
%        grayscale CT image volume (Analyze 7.5 format)
%        binary-valued segmentation image header (Analyze 7.5 format)
%        binary-valued segmentation image volume (Analyze 7.5 format)
%        log file with (x,y,z) coordinates of reference points (can be
%           generated by Analyze)
%        MR image header (Analyze 7.5 format)
%        MR image volume (Analyze 7.5 format)
%
%        text file with (x,y,z) coordinates of reference points
%        text file with azimuth, elevation, and radius of landmark points
%
%        standard input - not used
%
%     Outputs (located in "outputs_directory"):
%        .mat file with output variables from entire run
%        .csv file with output variables from entire run
%        .png files for (1) transverse section, (2) sagittal section, (3)
%           coronal section, and (4) surface normal attentuation profile.
%           The .png files are computed for each reference point.
%        standard output - not used
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
%     normalcy requires that the sampling of the input images in x, y, and
%     z be isotropic.
%
%     Code dependencies are indicated by the level of indent:
%     normalcy
%        calc_rot_angles
%           read_ref_landmark_file
%        rotate_3d_ct
%        calc_border_bbox
%        calc_landmarks
%           read_ref_landmarks_file
%        keep_connected_only
%        disp_transverse_movie
%        disp_sagittal_movie
%        disp_coronal_movie
%        seg_thresh
%        fminsearch (from Optimization Toolbox)
%           fit_function_rotation
%        rotate_3d_mr
%        fminsearch (from Optimization Toolbox)
%           fit_function_sphere
%        disp_transverse_normal
%        disp_sagittal_normal
%        disp_coronal_normal
%        find_edges_prof
%        find_edges_prof_mr
%        disp_transverse_normal_mr_csf
%        disp_sagittal_normal_mr_csf
%        disp_coronal_normal_mr_csf
%        Calc_braincase
%           calc_circumf_angles
%           rotate_ct_circumf
%        suptitle.m
%
%     normalcy is a function that is called by another function of Matlab
%     script or is invoked at the command line.
%
%  REVISION HISTORY
%     Version      Date                          Comment
%     -------   -----------      -------------------------------------------
%       1.0     12 May 2011       Initial release.
%       1.1     18 May 2011       Adding auto landmark calculation
%       1.2     23 June 2011      Improved image fill
%       1.3     7  July 2011      Added braincase segmentation 
%       1.4     21 July 2011      Added MR data and analyses
%       1.5     11 August 2011    Added 3D rotation and segmentation of
%                                 original grayscale CT image
%       1.6     1  September 2011 Added CSF measurement from MR data
%       1.7     8  September 2011 Added 3D rotation of original MR image
%                                 and alignment with CT image
%       1.8     27 October 2011   Added 3D registration of MR image with 
%                                 CT image
%       1.9     23 November 2011  Revised algorithm
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
%
function [] = normalcy(inputs_directory, outputs_directory, ...
    gray_file_name_root, ref_local_file_name, ...
    ref_landmark_file_name, out_file_name_root, ...
    mr_file_name_root, ...
    mr_ref_landmark_file_name, r_patch_mm, ...
    display_verbose, pause_before, pause_between)

    disp(' ') ;
    disp(' ') ;

    % Create the output directory.
    if (exist(outputs_directory, 'dir') == 0)
        [success, ~, ~] = mkdir(outputs_directory) ;
        if (success == 0)
            disp('ERROR: UNABLE TO CREATE OUTPUT DIRECTORY. EXITING!') ;
            return ;
        end
    end

    % Concatenate the directory and file names to find the complete path.
    
    full_ref_local_file_name = [inputs_directory '\' ...
        ref_local_file_name] ;
    full_ref_landmark_file_name = [inputs_directory '\' ...
        ref_landmark_file_name] ;
    full_gray_file_hdr_name  = ...
        [inputs_directory '\' gray_file_name_root  '.hdr'] ;
    full_gray_file_img_name = ...
        [inputs_directory '\' gray_file_name_root '.img'] ;
    full_mr_file_img_name = ...
        [inputs_directory '\' mr_file_name_root '.img'] ;
    full_mr_file_hdr_name  = ...
        [inputs_directory '\' mr_file_name_root  '.hdr'] ;
    full_mr_ref_landmark_file_name = [inputs_directory '\' ...
        mr_ref_landmark_file_name] ;

    % check to see if an MR file is present and set flag to 0 if present
    mr_empty = isempty(mr_file_name_root) ;
    
    % Read the MR image
    if (mr_empty == 0)
        mr_header = analyze75info(full_mr_file_hdr_name) ;
        mr_drcp_mm = [double(mr_header.PixelDimensions(1)) ...
                   double(mr_header.PixelDimensions(2)) ...
                   double(mr_header.PixelDimensions(3))] ;
        %mr_r_patch_rcp = r_patch_mm / mr_drcp_mm(2) ;
        mr_3d = double(analyze75read(full_mr_file_img_name)) ;
        
    % Calculate a threshold value for MR image. Find the mean intensity for
    % all voxels where intensity greater than 20 (assumed background) and
    % then divide by 2 to account for partial volume averaging.
    mr_head_thres = (mean(mr_3d(find(mr_3d > 20)))/2);    
    end
        
    % Read the grayscale and segmentation images. drcp_mm contains the voxel
    % dimensions. Assumed to be isotropic so only the Y dimension is read.
    % Adjusting to account for auto segmentation
    disp('Reading headers and image data ...') ;
    header = analyze75info(full_gray_file_hdr_name) ;
    drcp_mm = [double(header.PixelDimensions(1)) ...
               double(header.PixelDimensions(2)) ...
               double(header.PixelDimensions(3))] ;
    r_patch_rcp = r_patch_mm / drcp_mm(2) ;
    
    % Read in grayscale CT image
    g_3d = double(analyze75read(full_gray_file_img_name)) ;

   % Rotate grayscale CT image with its PAR as the reference. Return s_3d
   % which becomes the aligned and clipped binary file; g_3d which is the
   % aligned grayscale file; g_3d_orig is the original grayscale file
   % before alignment. Calc_rot_angles becomes the inputs for rotate_3d_ct.
   % Rotate_3d_ct does a connect operation so stray voxels should be gone,
   % but it does not do a fill operation.
    [image_center_in, alpha, beta, gamma, nas, par, pal]= ...
        calc_rot_angles(full_ref_local_file_name);
    image_center_out = image_center_in ;
    [g_3d_orig, s_3d, g_3d]= rotate_3d_ct(g_3d, image_center_in, image_center_out, ...
        alpha, beta, gamma);
    
    %Save aligned CT image
    full_align_file_name = [outputs_directory '\' out_file_name_root ...
        '.alignCT.vol'] ;
    fid=fopen(full_align_file_name,'w');
    g_3d_an75 = permute(g_3d,[2 1 3]);
    fwrite(fid,g_3d_an75,'int16');
    fclose(fid); 
    
    % set to even voxel value for later and saving
    par = round(par); 
    pal=round(pal);
    nas=round(nas);
    disp(['Aligned CT PAR = ' num2str(par)]) ;
    disp(['Aligned CT PAL = ' num2str(pal)]) ; 
    disp(['Aligned CT NAS = ' num2str(nas)]) ;

    % For CT data
    % Set the first two and last two rows, columns, and planes of the
    % segmentation image to zero. Some of the algorithms used below require
    % 2 "guard planes" on all sides.
    disp( ...
        'Creating 2 "guard" planes on all sides of the segmentation ...') ;
    s_3d(1:2,       :,         :        ) = 0 ;
    s_3d(end-1:end, :,         :        ) = 0 ;
    s_3d(:,         1:2,       :        ) = 0 ;
    s_3d(:,         end-1:end, :        ) = 0 ;
    s_3d(:,         :,         1:2      ) = 0 ;
    s_3d(:,         :,         end-1:end) = 0 ;
    
    % Clean up the segmentation image by filling in holes. This is
    % intentionally done in 2D in transverse, sagittal, and coronal planes. 
    temp0 = size(s_3d) ;
    num_row = temp0(1) ;
    num_col = temp0(2) ;
    num_pln = temp0(3) ;
    
    disp(['Filling in holes of segmentation image on 2D transverse ' ...
        'slices ...']) ;
    temp1 = zeros(size(s_3d)) ;
    parfor pln_num = 1:num_pln
        temp1(:, :, pln_num) = imfill(s_3d(:, :, pln_num), 'holes') ;
    end
    
    disp(['Filling in holes of segmentation image on 2D sagittal ' ...
        'slices ...']) ;
    temp2 = zeros(size(temp1)) ;
    parfor col_num = 1:num_col
        temp2(:, col_num, :) = imfill(squeeze(temp1(:, col_num, :)), 'holes') ;
    end
    
    disp(['Filling in holes of segmentation image on 2D coronal ' ...
        'slices ...']) ;
    temp3 = zeros(size(temp2)) ;
    parfor row_num = 1:num_row
        temp3(row_num, :, :) = imfill(squeeze(temp2(row_num, :, :)), 'holes') ;
    end
    
    s_3d = temp3 ;
    
    % Calculate a bounding box.   
    [min_rcp_of_data, max_rcp_of_data, list_border_bc_rcp, b_3d_bc] = ...
            calc_border_bbox(s_3d);
    % At this point list_border_bc_rcp is the border points or the segmented
    % and filled CT head. 

% Only if MR data
if (mr_empty == 0)
    
    % CT registation landmark calculations only -------------------------------------

    % Auto calulate landmarks based on spherical coordinate system. This
    % function reads in the local coordinate system points, reads in the
    % specified azimuth and elevation angles, calculates a line and border
    % voxel intersection, and returns the inputs to the main loop below.
    % Border voxels must have been calculated prior to this function call.

       [ref_point_mat_xyz, num_ref_points, local_origin] = ...
        calc_landmarks(full_mr_ref_landmark_file_name, list_border_bc_rcp, ...
                       min_rcp_of_data, max_rcp_of_data, nas, par, pal) ;
                                
        for nn=1:num_ref_points
            disp(['CT ref pt = ' num2str(ref_point_mat_xyz(nn,:))]) ;
        end 
%--------------------------------------------------------------------        
%         %Display Registration landmarks on CT image
%         hold on; % for figure(1)
%         ref_point_mat_rcp = [ref_point_mat_xyz(:, 2) ...
%             ref_point_mat_xyz(:, 1) ref_point_mat_xyz(:, 3)];
%         ref_pt_empty = isempty(ref_point_mat_xyz) ;
%         if (ref_pt_empty == 0)
%             temp1 = size(ref_point_mat_rcp) ;
%             land_tot = temp1(1) ;
%             parfor land = 1:land_tot
%                 hold on;
%                 coord3d = ref_point_mat_rcp(land,:) ;
%                 xc = coord3d(1) ;
%                 yc = coord3d(2) ;
%                 zc = coord3d(3) ;
%                 [x,y,z]=ellipsoid(xc,yc,zc,7,7,7) ;
%                 mesh(x,y,z) ;
%             end
%         else
%             disp('ERROR:  NO REFERENCE POINTS SPECIFIED. EXITING!') ;
%             return ;   
%         end
        
   %------------------------------------------------------------------
   % For Registration calculate the landmarks and patches of voxels
    
    % Process all of the reference points. This is the main loop of this
    % function. This will be for the registration points only.
    
    % num_rand determines the number of points to be sampled within
    % each patch, therefore, the total number of ct ref pts is rand_num
    % times the number of landmarks used, so best to keep small for
    % runtime issues.
    % Does fewer random pts work better? Should the number of pts be
    % increased per iteration?
   
    num_rand = 10 ;  % Defines number of points to randomly sample from a patch
%    cum_ct_pts_rcp = zeros((num_rand * num_ref_points), 3);
    cum_ct_pts_rcp = [];
    num_reg_pts = num_rand * num_ref_points ;
    % put parfor here
    parfor ref_point_numb = 1:num_ref_points
        ref_point_xyz = ref_point_mat_xyz(ref_point_numb, :) ;
        ref_point_rcp = [ref_point_xyz(2) ref_point_xyz(1) ...
                         ref_point_xyz(3)] ;
        disp(' ') ;
        disp(['Processing reference point ' num2str(ref_point_numb) ...
            ' out of ' num2str(num_ref_points) '.']) ;
        
        % Clean up the border by discarding all stray voxels that are not
        % connected to the actual border. First, find a "seed point," a
        % point on the borderm which is closest to the reference point.
        
%       It might be possible to delete this next section but need to check
%       the impact of all the variables like b_3d_bc. 
        disp('Cleaning up border by discarding unconnected voxels ...') ;
        
        dist_away_rcp = ...
            sqrt(((list_border_bc_rcp(:,1) - ref_point_rcp(1)).^2) ...
               + ((list_border_bc_rcp(:,2) - ref_point_rcp(2)).^2) ...
               + ((list_border_bc_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        [~, index_of_closest_point] = min(dist_away_rcp) ;
        rcp_seed = [list_border_bc_rcp(index_of_closest_point, 1) ...
                    list_border_bc_rcp(index_of_closest_point, 2) ...
                    list_border_bc_rcp(index_of_closest_point, 3)] ;
        
        [b_3d, list_border_rcp] = keep_connected_only(b_3d_bc, rcp_seed) ;
        disp('*********************************************************') ;
        disp(['Number of border voxels before cleanup = ' ...
            num2str(sum(sum(sum(b_3d_bc)))) '.']) ;
        disp(['Number of border voxels after cleanup =  ' ...
            num2str(sum(sum(sum(b_3d)))) '.']) ;
        disp('*********************************************************') ;
%
        % Now define the "patch" of points on the border that will be used
        % to fit the sphere. First, we eliminate all points that are too
        % far away from the reference point.

        disp(['Creating a "patch" of voxels that the sphere will be ' ...
              'fitted to ...']) ;
        dist_away_rcp = sqrt( ...
              ((list_border_rcp(:,1) - ref_point_rcp(1)).^2) ...
            + ((list_border_rcp(:,2) - ref_point_rcp(2)).^2) ...
            + ((list_border_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        rrr = find(dist_away_rcp <= r_patch_rcp) ;
        list_patch_bc_rcp = [] ;
        list_patch_bc_rcp(:,1) = list_border_rcp(rrr, 1) ;
        list_patch_bc_rcp(:,2) = list_border_rcp(rrr, 2) ;
        list_patch_bc_rcp(:,3) = list_border_rcp(rrr, 3) ;

        % Build an image volume for the patch. At this time, the patch is
        % the set of all border voxels which are close enough to the
        % reference point.

        p_3d_bc = zeros(size(s_3d)) ;
        for i = 1:length(list_patch_bc_rcp(:,1))
            p_3d_bc(list_patch_bc_rcp(i,1), ...
                    list_patch_bc_rcp(i,2), ...
                    list_patch_bc_rcp(i,3)) = 1 ;
        end
        
        % Clean up the patch by discarding all stray voxels that are not
        % connected to the seeds point. The seed point is the voxel on the
        % border closest to the reference point.

        disp('Cleaning up patch by discarding unconnected voxels ...') ;
        dist_away_rcp = ...
            sqrt(((list_patch_bc_rcp(:,1) - ref_point_rcp(1)).^2) ...
               + ((list_patch_bc_rcp(:,2) - ref_point_rcp(2)).^2) ...
               + ((list_patch_bc_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        [~, index_of_closest_point] = min(dist_away_rcp) ;

        rcp_seed = [list_patch_bc_rcp(index_of_closest_point, 1) ...
                    list_patch_bc_rcp(index_of_closest_point, 2) ...
                    list_patch_bc_rcp(index_of_closest_point, 3)] ;

        [p_3d, list_patch_rcp] = keep_connected_only(p_3d_bc, rcp_seed) ;
 
        patch_size = (sum(sum(sum(p_3d))));
        rand_sel = randi(patch_size, 1, num_rand) ; 
        ct_pts_rcp = zeros(num_rand, 3) ;
        for cntr=1:num_rand
            ct_pts_rcp(cntr, :) = list_patch_rcp(rand_sel(cntr), :);
        end

        cum_ct_pts_rcp = [cum_ct_pts_rcp; ct_pts_rcp] ;
        
        disp('*********************************************************') ;
        disp(['Number of patch voxels before cleanup = ' ...
            num2str(sum(sum(sum(p_3d_bc)))) '.']) ;
        disp(['Number of patch voxels after cleanup =  ' ...
            num2str(sum(sum(sum(p_3d)))) '.']) ;
        % this needs fixed to display the voxels kept
        disp(['Number of voxels kept for registration =  ' ...
            num2str(num_rand) '.']) ;
        disp('*********************************************************') ;

        % Build an image volume that is equal to one at the border voxels
        % and equal to two at the border voxels that are close enough to
        % the reference point. This image is used for display purposes
        % only.

        b_3d_with_patch = b_3d + p_3d ;

        % Display all of the planes of the images.

        if (display_verbose == 1)
            disp('Displaying a movie of transverse planes ...') ;
            disp_transverse_movie(g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
            disp('Displaying a movie of sagittal planes ...') ;
            disp_sagittal_movie(  g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
            disp('Displaying a movie of coronal planes ...') ;
            disp_coronal_movie(   g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
        end
    end
    
    % This is the end of the loop that calculates the patches for
    % the ct registration points. Each patch is randomly sampled and 
    % the cummulative points are used as the ct refence points for the 
    % registraiton. The ref_points need to be
    % passed into the registration routine along with the mr boundary
    % points. The interstection between the ct line and the mr boundary
    % point becomes the mr ref point and is used to minimize the error.
    %------------------------------------------------------------------
end
    
    % Use this when not prealigning. Segment MR head, make binary, Display
    clip_plane=par(3);
    if (mr_empty == 0)
        [mr_s_3d] = seg_thresh(mr_3d, mr_head_thres, clip_plane) ;
               
        % For MR data
        % Set the first two and last two rows, columns, and planes of the
        % segmentation image to zero. Some of the algorithms used below require
        % 2 "guard planes" on all sides.
        disp('Creating 2 "guard" planes on all sides of the segmentation ...') ;
        mr_s_3d(1:2,       :,         :        ) = 0 ;
        mr_s_3d(end-1:end, :,         :        ) = 0 ;
        mr_s_3d(:,         1:2,       :        ) = 0 ;
        mr_s_3d(:,         end-1:end, :        ) = 0 ;
        mr_s_3d(:,         :,         1:2      ) = 0 ;
        mr_s_3d(:,         :,         end-1:end) = 0 ;
    
        % Clean up the segmentation image by filling in holes. This is
        % intentionally done in 2D in transverse, sagittal, and coronal planes.
        temp0 = size(mr_s_3d) ;
        num_row = temp0(1) ;
        num_col = temp0(2) ;
        num_pln = temp0(3) ;
    
        disp(['Filling in holes of segmentation image oncalc 2D transverse ' ...
            'slices ...']) ;
        temp1 = zeros(size(mr_s_3d)) ;
        parfor pln_num = 1:num_pln
            temp1(:, :, pln_num) = imfill(mr_s_3d(:, :, pln_num), 'holes') ;
        end
    
        disp(['Filling in holes of segmentation image on 2D sagittal ' ...
            'slices ...']) ;
        temp2 = zeros(size(temp1)) ;
        parfor col_num = 1:num_col
            temp2(:, col_num, :) = imfill(squeeze(temp1(:, col_num, :)), 'holes') ;
        end
    
        disp(['Filling in holes of segmentation image on 2D coronal ' ...
        'slices ...']) ;
        temp3 = zeros(size(temp2)) ;
        parfor row_num = 1:num_row
            temp3(row_num, :, :) = imfill(squeeze(temp2(row_num, :, :)), 'holes') ;
        end
    
        mr_s_3d = temp3 ;
   
        % Calculate a bounding box for the segmentation.
        
        [mr_min_rcp_of_data, mr_max_rcp_of_data, mr_list_border_bc_rcp, mr_b_3d_bc] = ...
            calc_border_bbox(mr_s_3d);

    end
 %-----------------------------------------------------------------
%
tic
    disp(['Beginning MR to CT registration']) ;
    % CT and MR registration, only if MR --------------------------------------
    if (mr_empty == 0) 

        % For rotation, initialize center to actual volume center and 0 deg
        % for rotation angles
        temp0 = size(mr_s_3d) ;
        num_x = temp0(2) ;
        num_y = temp0(1) ;
        num_z = temp0(3) ;
        
        % center_init_xyz and therefore center_xyz should be the centroid
        % of the binary mr data to speed up registration. probably not a
        % big deal for heads as the center of the volume is well within the
        % head and near the center of rotation anyway.

        step_size=460; %step is actually 5% of 1/2 this value
        center_xyz = [(num_x+1)/2 (num_y+1)/2 (num_z+1)/2]; % center of volume
        alpha_est = 0;
        beta_est = 0 ;
        gamma_est = 0 ;
        center_est_xyz(1) = 0 ;
        center_est_xyz(2) = 0 ;
        center_est_xyz(3) = 0 ;

        sse_prev = 0;
        reg_loop = 4;
        
        for iter = 1:reg_loop
            step_size = step_size/2 ;
                
             parms_init(1) = alpha_est + step_size ;
             parms_init(2) = beta_est + step_size ;
             parms_init(3) = gamma_est + step_size ;
             parms_init(4) = center_est_xyz(1) + step_size ;
             parms_init(5) = center_est_xyz(2) + step_size ;
             parms_init(6) = center_est_xyz(3) + step_size ;
        
            disp(['Starting registration iteration = ' num2str(iter)]) ;
         
            % pass in a random sampling of ct pts that come from the patches
            cum_ct_pts_xyz = cum_ct_pts_rcp(:,[2 1 3]);
           
            % Estimate the parameters of the best-fit rotation and translation
            % using fminsearch, passing in MR border points. 
            %Will calculate corresponding points within the function.
            % It would be best if we clipped off the bottom portion of
            % mr_list_border_bc_rcp prior to sending it in. Could do this
            % with lowest z value of bounding box?
            f = @(x)fit_function_rotation(x, cum_ct_pts_xyz, ...
                        mr_list_border_bc_rcp, center_xyz, ...
                        local_origin, step_size) ;
            [parms_out, sse_reg, exitflag, foutput] = fminsearch(f, parms_init, optimset('MaxFunEvals',4000)) ;
            
             disp(['Sum of squared errors between Registered landmarks = ' ...
                num2str(sse_reg)]) ;
            if (exitflag == 1)
                disp(['Optimization routine reported satisfactory ' ...
                      'convergence.']) ;
            else
                disp(['WARNING: OPTIMIZATION ROUTINE FAILED. EXITFLAG = ' ...
                    num2str(exitflag) '.']) ;
            end
               
           % Adjust and save parameters
            center_est_xyz(1:3) = parms_out(4:6) - step_size;
            alpha_est      = parms_out(1) - step_size ;
            beta_est      = parms_out(2) - step_size ;
            gamma_est      = parms_out(3) - step_size ;
            
            if (or(abs(sse_reg - sse_prev) < 1, iter == reg_loop))
                break;
            end
            sse_prev = sse_reg;   
        end
        
            %Rotate MR based on fit results center_init_xyz or center_xyz
            center_est_xyz(1:3) = parms_out(4:6) + center_xyz - step_size;
            clip_plane = par(3);  %CT PAR plane
            [mr_3d_orig, mr_s_3d, mr_3d]= rotate_3d_mr(mr_3d, center_xyz, center_est_xyz, ...
                alpha_est, beta_est, gamma_est, clip_plane, mr_head_thres);   
            
            % For MR data
            % Set the first two and last two rows, columns, and planes of the
            % segmentation image to zero. Some of the algorithms used below require
            % 2 "guard planes" on all sides.
            disp('Creating 2 "guard" planes on all sides of the segmentation ...') ;
            mr_s_3d(1:2,       :,         :        ) = 0 ;
            mr_s_3d(end-1:end, :,         :        ) = 0 ;
            mr_s_3d(:,         1:2,       :        ) = 0 ;
            mr_s_3d(:,         end-1:end, :        ) = 0 ;
            mr_s_3d(:,         :,         1:2      ) = 0 ;
            mr_s_3d(:,         :,         end-1:end) = 0 ;

            % Clean up the segmentation image by filling in holes. This is
            % intentionally done in 2D in transverse, sagittal, and coronal planes.
            temp0 = size(mr_s_3d) ;
            num_row = temp0(1) ;
            num_col = temp0(2) ;
            num_pln = temp0(3) ;

            disp(['Filling in holes of segmentation image on 2D transverse ' ...
                'slices ...']) ;
            temp1 = zeros(size(mr_s_3d)) ;
            parfor pln_num = 1:num_pln
                temp1(:, :, pln_num) = imfill(mr_s_3d(:, :, pln_num), 'holes') ;
            end

            disp(['Filling in holes of segmentation image on 2D sagittal ' ...
                'slices ...']) ;
            temp2 = zeros(size(temp1)) ;
            parfor col_num = 1:num_col
                temp2(:, col_num, :) = imfill(squeeze(temp1(:, col_num, :)), 'holes') ;
            end

            disp(['Filling in holes of segmentation image on 2D coronal ' ...
            'slices ...']) ;
            temp3 = zeros(size(temp2)) ;
            parfor row_num = 1:num_row
                temp3(row_num, :, :) = imfill(squeeze(temp2(row_num, :, :)), 'holes') ;
            end

            mr_s_3d = temp3 ;

            % need to recalculate border points and bounding box after rotation
            [mr_min_rcp_of_data, mr_max_rcp_of_data, mr_list_border_bc_rcp, mr_b_3d_bc] = ...
                calc_border_bbox(mr_s_3d);
            
        % MR landmark calculations. Only used for visualization.
        % Auto calulate landmarks based on spherical coordinate system. This
        % function reads in the local coordinate system points, reads in the
        % specified azimuth and elevation angles, calculates a line and border
        % voxel intersection, and returns the inputs to the main loop below.
        % Border voxels must have been calculated prior to this function call.
        % Don't clip MR volume until after registered otherwise will lose some
        % of the volume during registration. The local coordinate system is
        % based on the CT data as that is what is used for profile
        % plotting.
           [mr_ref_point_mat_xyz, mr_num_ref_points, mr_local_z] = ...
            calc_landmarks(full_mr_ref_landmark_file_name, mr_list_border_bc_rcp, ...
            min_rcp_of_data, max_rcp_of_data, nas, par, pal) ;
            for nn=1:mr_num_ref_points
                disp(['MR ref pt = ' num2str(mr_ref_point_mat_xyz(nn,:))]) ;
            end

         % Plotting position of landmarks on segmented displayed surface.
         % The points are plotted in rcp to account for some anomaly in
         % displaying on top of the volume image. This could be due to axes
         % handles. The ellipsoid/mesh plot axes are likely different from
         % the default figure axes and inherit the wrong axes.
%          figure(1); %to display on CT head, remove if wanted on new MR data
%             hold on;
%             mr_ref_point_mat_rcp = [mr_ref_point_mat_xyz(:, 2) ...
%             mr_ref_point_mat_xyz(:, 1) mr_ref_point_mat_xyz(:, 3)];
%             ref_pt_empty = isempty(mr_ref_point_mat_xyz) ;
%             if (ref_pt_empty == 0)
%                 temp1 = size(mr_ref_point_mat_rcp) ;
%                 land_tot = temp1(1) ;
%                     parfor land = 1:land_tot
%                         hold on;
%                         coord3d = mr_ref_point_mat_rcp(land,:) ;
%                      xc = coord3d(1) ;
%                      yc = coord3d(2) ;
%                         zc = coord3d(3) ;
%                         [x,y,z]=ellipsoid(xc,yc,zc,7,7,7) ;
%                         mesh(x,y,z) ;
%                     end
%             else
%                disp('ERROR:  NO REFERENCE POINTS SPECIFIED. EXITING!') ; 
%                return;
%             end  

    %Save aligned MR image
        full_mr_align_file_name = [outputs_directory '\' out_file_name_root ...
            '.alignMR.vol'] ;
        fid=fopen(full_mr_align_file_name,'w');
        mr_3d_an75 = permute(mr_3d,[2 1 3]);
        fwrite(fid,mr_3d_an75,'int16');
        fclose(fid);
        
   end
toc
 %-------------------------------------------------------------------------       
    % Auto calulate landmarks based on spherical coordinate system. This
    % function reads in the local coordinate system points, reads in the
    % specified azimuth and elevation angles, calculates a line and border
    % voxel intersection, and returns the inputs to the main loop below.
    % Border voxels must have been calculated prior to this function call.
       
    [ref_point_mat_xyz, num_ref_points, local_origin] = ...
        calc_landmarks(full_ref_landmark_file_name, list_border_bc_rcp, ...
                       min_rcp_of_data, max_rcp_of_data, nas, par, pal) ;       
    
    %____Display Only___________________________________________________
    %Display reference landmarks on CT image
%         hold on; % for figure(1)
%         ref_point_mat_rcp = [ref_point_mat_xyz(:, 2) ...
%             ref_point_mat_xyz(:, 1) ref_point_mat_xyz(:, 3)];
%         ref_pt_empty = isempty(ref_point_mat_xyz) ;
%         if (ref_pt_empty == 0)
%             temp1 = size(ref_point_mat_rcp) ;
%             land_tot = temp1(1) ;
%             parfor land = 1:land_tot
%                 hold on;
%                 coord3d = ref_point_mat_rcp(land,:) ;
%                 xc = coord3d(1) ;
%                 yc = coord3d(2) ;
%                 zc = coord3d(3) ;
%                 [x,y,z]=ellipsoid(xc,yc,zc,7,7,7) ;
%                 mesh(x,y,z) ;
%             end
%         else
%             disp('ERROR:  NO REFERENCE POINTS SPECIFIED. EXITING!') ;
%             return ;   
%         end
    %___________________________________________________________________
    % Initialize the vectors that will be estimated.

    center_est_mat_xyz    = zeros(num_ref_points, 3) ;
    radius_est_vec_mm     = zeros(num_ref_points, 1) ;
    direction_vec_mat_xyz = zeros(num_ref_points, 3) ;
    
 
    % Process all of the reference points. This is the main loop of this
    % function.
    
    for ref_point_num = 1:num_ref_points
        ref_point_xyz = ref_point_mat_xyz(ref_point_num, :) ;
        ref_point_rcp = [ref_point_xyz(2) ref_point_xyz(1) ...
                         ref_point_xyz(3)] ;
        disp(' ') ;
        disp(['Processing reference point ' num2str(ref_point_num) ...
            ' out of ' num2str(num_ref_points) '.']) ;
        
        % Clean up the border by discarding all stray voxels that are not
        % connected to the actual border. First, find a "seed point," a
        % point on the borderm which is closest to the reference point.

        disp('Cleaning up border by discarding unconnected voxels ...') ;
        dist_away_rcp = ...
            sqrt(((list_border_bc_rcp(:,1) - ref_point_rcp(1)).^2) ...
               + ((list_border_bc_rcp(:,2) - ref_point_rcp(2)).^2) ...
               + ((list_border_bc_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        [~, index_of_closest_point] = min(dist_away_rcp) ;
        rcp_seed = [list_border_bc_rcp(index_of_closest_point, 1) ...
                    list_border_bc_rcp(index_of_closest_point, 2) ...
                    list_border_bc_rcp(index_of_closest_point, 3)] ;
        [b_3d, list_border_rcp] = keep_connected_only(b_3d_bc, rcp_seed) ;
        disp('*********************************************************') ;
        disp(['Number of border voxels before cleanup = ' ...
            num2str(sum(sum(sum(b_3d_bc)))) '.']) ;
        disp(['Number of border voxels after cleanup =  ' ...
            num2str(sum(sum(sum(b_3d)))) '.']) ;
        disp('*********************************************************') ;

        % Now define the "patch" of points on the border that will be used
        % to fit the sphere. First, we eliminate all points that are too
        % far away from the reference point.

        disp(['Creating a "patch" of voxels that the sphere will be ' ...
              'fitted to ...']) ;
        dist_away_rcp = sqrt( ...
              ((list_border_rcp(:,1) - ref_point_rcp(1)).^2) ...
            + ((list_border_rcp(:,2) - ref_point_rcp(2)).^2) ...
            + ((list_border_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        rrr = find(dist_away_rcp <= r_patch_rcp) ;
        list_patch_bc_rcp = [] ;
        list_patch_bc_rcp(:,1) = list_border_rcp(rrr, 1) ;
        list_patch_bc_rcp(:,2) = list_border_rcp(rrr, 2) ;
        list_patch_bc_rcp(:,3) = list_border_rcp(rrr, 3) ;

        % Build an image volume for the patch. At this time, the patch is
        % the set of all border voxels which are close enough to the
        % reference point.

        p_3d_bc = zeros(size(s_3d)) ;
        for i = 1:length(list_patch_bc_rcp(:,1))
            p_3d_bc(list_patch_bc_rcp(i,1), ...
                    list_patch_bc_rcp(i,2), ...
                    list_patch_bc_rcp(i,3)) = 1 ;
        end
        
        % Clean up the patch by discarding all stray voxels that are not
        % connected to the seeds point. The seed point is the voxel on the
        % border closest to the reference point.

        disp('Cleaning up patch by discarding unconnected voxels ...') ;
        dist_away_rcp = ...
            sqrt(((list_patch_bc_rcp(:,1) - ref_point_rcp(1)).^2) ...
               + ((list_patch_bc_rcp(:,2) - ref_point_rcp(2)).^2) ...
               + ((list_patch_bc_rcp(:,3) - ref_point_rcp(3)).^2)) ;
        [~, index_of_closest_point] = min(dist_away_rcp) ;

        rcp_seed = [list_patch_bc_rcp(index_of_closest_point, 1) ...
                    list_patch_bc_rcp(index_of_closest_point, 2) ...
                    list_patch_bc_rcp(index_of_closest_point, 3)] ;

        [p_3d, list_patch_rcp] = keep_connected_only(p_3d_bc, rcp_seed) ;

        list_patch_rcp(:,1) = list_patch_rcp(:,1) ;
        list_patch_rcp(:,2) = list_patch_rcp(:,2) ;
        list_patch_rcp(:,3) = list_patch_rcp(:,3) ;
        disp('*********************************************************') ;
        disp(['Number of patch voxels before cleanup = ' ...
            num2str(sum(sum(sum(p_3d_bc)))) '.']) ;
        disp(['Number of patch voxels after cleanup =  ' ...
            num2str(sum(sum(sum(p_3d)))) '.']) ;
        disp('*********************************************************') ;

        % Build an image volume that is equal to one at the border voxels
        % and equal to two at the border voxels that are close enough to
        % the reference point. This image is used for display purposes
        % only.

        b_3d_with_patch = b_3d + p_3d ;

        % Display all of the planes of the images.

        if (display_verbose == 1)
            disp('Displaying a movie of transverse planes ...') ;
            disp_transverse_movie(g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
            disp('Displaying a movie of sagittal planes ...') ;
            disp_sagittal_movie(  g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
            disp('Displaying a movie of coronal planes ...') ;
            disp_coronal_movie(   g_3d, s_3d, b_3d_with_patch, ...
                pause_before, pause_between) ;
        end

        % Estimate the parameters of the best-fit sphere using fminsearch.
        
        disp('Performing fminsearch to find sphere of best fit ...') ;
        center_init_rcp(1) = 0.5 * ...
            (min_rcp_of_data(1) + max_rcp_of_data(1)) ;
        center_init_rcp(2) = 0.5 * ...
            (min_rcp_of_data(2) + max_rcp_of_data(2)) ;
        center_init_rcp(3) = 0.5 * ...
            (min_rcp_of_data(3) + max_rcp_of_data(3)) ;
        radius_init_rcp = mean(sqrt( ...
              ((list_patch_rcp(:,1) - center_init_rcp(1)).^2) ...
            + ((list_patch_rcp(:,2) - center_init_rcp(2)).^2) ...
            + ((list_patch_rcp(:,3) - center_init_rcp(3)).^2))) ;
        parms_init = [center_init_rcp(1) center_init_rcp(2) ...
                      center_init_rcp(3) radius_init_rcp] ;

        f = @(x)fit_function_sphere(x, list_patch_rcp(:,1), ...
            list_patch_rcp(:,2), list_patch_rcp(:,3)) ;
        [parms_out, sse_search, exitflag] = fminsearch(f, parms_init) ;

        center_est_rcp(1:3) = parms_out(1:3) ;
        radius_est_rcp      = parms_out(4)   ;

        center_est_xyz      = center_est_rcp([2 1 3]) ;
        center_est_mat_xyz(ref_point_num, :) = center_est_xyz ;

        radius_est_mm       = radius_est_rcp * drcp_mm(2) ;
        radius_est_vec_mm(ref_point_num) = radius_est_mm ;

        % Now display the results if CT only.
        if (mr_empty ~= 0)
            disp(['Displaying the transverse plane at the reference point ' ...
                'and the normal vector ...']) ;
            disp_transverse_normal(g_3d, b_3d, p_3d, ref_point_rcp, ...
                center_est_rcp, pause_before) ;
            figure(4) ;
            title_string = ['Reference Point xyz = (' ...
                num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ...
                ',' num2str(ref_point_xyz(3)) ')'] ;
            suptitle(title_string) ;
            print_string = ['print -dpng -r600 ' outputs_directory '\' ...
                out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
                '_transverse.png'] ;
            eval(print_string) ;
            disp(['Displaying the sagittal plane at the reference point ' ...
                'and the normal vector ...']) ;
            disp_sagittal_normal(g_3d, b_3d, p_3d, ref_point_rcp, ...
                center_est_rcp, pause_before) ;
            figure(5) ;
            title_string = ['Reference Point xyz = (' ...
                num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
                num2str(ref_point_xyz(3)) ')'] ;
            suptitle(title_string) ;
            print_string = ['print -dpng -r600 ' outputs_directory '\' ...
                out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
                '_sagittal.png'] ;
            eval(print_string) ;
            disp(['Displaying the coronal plane at the reference point ' ...
                'and the normal vector ...']) ;
            disp_coronal_normal(g_3d, b_3d, p_3d, ref_point_rcp, ...
                center_est_rcp, pause_before) ;
            figure(6) ;
            title_string = ['Reference Point xyz = (' ...
                num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
                num2str(ref_point_xyz(3)) ')'] ;
            suptitle(title_string) ;
            print_string = ['print -dpng -r600 ' outputs_directory '\' ...
                out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
                '_coronal.png'] ;
            eval(print_string) ;
        end
        
        % Now find the profile of the grayscale image along the line from
        % the reference point to the center of the circle that fit it best.

        disp('Computing the direction vector for the normal vector ...') ;
        temp = [(center_est_rcp(1) - ref_point_rcp(1))   ...
                (center_est_rcp(2) - ref_point_rcp(2))   ...
                (center_est_rcp(3) - ref_point_rcp(3))] ;
        dist_ref_to_cent = norm(temp) ;
        direction_vec_rcp = temp / dist_ref_to_cent ;
        direction_vec_xyz    = direction_vec_rcp([2 1 3]) ;
        direction_vec_mat_xyz(ref_point_num, :) = direction_vec_xyz ;

        disp(['Displaying the profile through the CT image along the ' ...
            'normal vector ...']) ;
        axis_x_min_mm = -10 ;
        axis_x_max_mm =  30 ;
       
        delta_interp_mm = 0.01 ;
        range_mm = axis_x_min_mm:delta_interp_mm:axis_x_max_mm ;
        range_vox = range_mm / drcp_mm(2) ;

        if (ref_point_num == 1)
            prof_mat = zeros(num_ref_points, length(range_mm)) ;
        end

        row_temp = ref_point_rcp(1) + (range_vox * direction_vec_rcp(1)) ;
        col_temp = ref_point_rcp(2) + (range_vox * direction_vec_rcp(2)) ;
        pln_temp = ref_point_rcp(3) + (range_vox * direction_vec_rcp(3)) ;
        prof = interp3(g_3d, col_temp, row_temp, pln_temp) ;  % CT data
        kkk = find(isnan(prof)) ;
        prof(kkk) = -1000 ;
        prof_mat(ref_point_num, :) = prof ;
        
        % values to classify tissues at 0 crossings (peaks) according to tissue type
        tissue_type = 1:6 ;

        % set thresholds for tissues [skin, fat, muscle, trabecula, cortical].
        % Adjusted values of skin and trabecula to account for noise. Look at
        % adding more intelligence into algorithm in find edges.
        % tissue_thr = [-750 -200 -25 200 600] ;
        tissue_thr = [-550 -200 -25 200 600] ;
        skin_thr = tissue_thr(1) ;
        fat_thr = tissue_thr(2) ;
        mus_thr = tissue_thr(3) ;
        trab_thr = tissue_thr(4) ;
        cort_thr = tissue_thr(5) ;

        % Create an nx3 matrix to store values of ST_mm, BT_mm, and Ave_HU_Bone
        % where n = number of landmark points being examined. Initializing
        % the vector prior to calling find_edges_prof to ensure 0's in case
        % find_edges_prof fails.

        if (ref_point_num == 1)
           find_prof_results = zeros(num_ref_points, 3) ;
        end
        
        % Find values in profile plot where 1st derivitave crosses zero
        [prof_edges, prof_edges_new, prof_edges_new_index, first_pt, ...
            skin_out2_pt, bone_out1_pt, bone_out2_pt, bone_in1_pt, bone_in2_pt, ...
            sk_edge_i, bo_edge_i, bi_edge_i, Ave_HU_Bone, ST_mm, BT_mm, steps_error_vec] = ...
            find_edges_prof(prof, range_mm, tissue_type, tissue_thr) ;
        
        % Upon return from find_edges_prof check for errors. Current total = 17
        if sum(steps_error_vec) ~= 17 ;
            continue ;
        end

        % store values in a matrix indexed by num_ref_points. This allows the
        % values to be written to the csv file at the end of normalcy. 
        find_prof_results(ref_point_num, 1) = ST_mm ;
        find_prof_results(ref_point_num, 2) = BT_mm ;
        find_prof_results(ref_point_num, 3) = Ave_HU_Bone ;
        % store the bone_inner HU values for each landmark for use in braincase
        % volume segmentation
        bi_edge_HU_values(ref_point_num) = prof(bi_edge_i) ; 

        figure(7) ;
        % plot profile on top and gradients underneath
            subplot(1,1,1), plot(range_mm, prof) ;
            axis_y_min_mm =  -1025 ;
            axis_y_max_mm = (floor(max(prof) / 100) + 1) * 100 ;
            axis([axis_x_min_mm axis_x_max_mm axis_y_min_mm axis_y_max_mm]) ;
            xlabel('Distance from Reference Point Along Surface Normal (mm)') ;
            ylabel('Image Intensity (HU)') ;

        % Display HU levels on plots
            hold on ;
            plot([axis_x_min_mm axis_x_max_mm],[skin_thr skin_thr],'--','Color',[1 0 1]) ;
            hold on ; 
            plot([axis_x_min_mm axis_x_max_mm],[fat_thr fat_thr],'--yellow') ;
            hold on ;
            plot([axis_x_min_mm axis_x_max_mm],[mus_thr mus_thr],'--red') ;
            hold on ;
            plot([axis_x_min_mm axis_x_max_mm],[trab_thr trab_thr],'--green') ;
            hold on ;
            plot([axis_x_min_mm axis_x_max_mm],[cort_thr cort_thr],'--blue') ;             

        % plot profile peaks based on zero crossing of gradient                      
            hold on ;

        % calculate x coordinate at which to plot HU at a zero crossing (peak)
             ccc1= -10.01 + (0.01 * prof_edges_new_index) ;        
             plot (ccc1, prof_edges_new, '+yellow') ; 

        %plot first point in blue
            hold on ;
            plot (ccc1(first_pt), prof_edges_new(first_pt), '+blue') ;

        %plot second skin point in blue
            hold on ;
            plot (ccc1(skin_out2_pt), prof_edges_new(skin_out2_pt), '+blue') ;        

        %plot first outer bone point in cyan
            hold on ;
            plot (ccc1(bone_out1_pt), prof_edges_new(bone_out1_pt), '+red') ;

        %plot second outer bone point in cyan
            hold on ;
            plot (ccc1(bone_out2_pt), prof_edges_new(bone_out2_pt), '+red') ;        

        %plot first inner bone point in green
            hold on ;
            plot (ccc1(bone_in1_pt), prof_edges_new(bone_in1_pt), '+green') ; 

        %plot second inner bone point in green
            hold on ;
            plot (ccc1(bone_in2_pt), prof_edges_new(bone_in2_pt), '+green') ;        

        %Plot a vertical line denoting the skin edge
            hold on ;
            plot ([range_mm(sk_edge_i) range_mm(sk_edge_i)], ...
                [axis_y_min_mm axis_y_max_mm], 'blue') ;

        %Plot a vertical line denoting the outer bone edge
            hold on ;
            plot ([range_mm(bo_edge_i) range_mm(bo_edge_i)], ...
                [axis_y_min_mm axis_y_max_mm], 'red') ;

        %Plot a vertical line denoting the inner bone edge
            hold on ;
            plot ([range_mm(bi_edge_i) range_mm(bi_edge_i)], ...
                [axis_y_min_mm axis_y_max_mm], 'green') ;       

        subplot(1,1,1) ;
        title_string = ['Subject ' num2str(gray_file_name_root(1:6)) ',' ...
            '  Reference Point xyz = (' ...
            num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
            num2str(ref_point_xyz(3)) ')'] ;
        title_string1 = ['Average Bone HU = ' num2str(Ave_HU_Bone)] ;
        title_string2 = ['Soft Tissue Thickness (mm) = ' num2str(ST_mm)...
            '  Bone Thickness (mm) = ' num2str(BT_mm)] ;
        title({title_string, title_string1, title_string2}) ;

        % This plots the Main Title at the very top. Comes after all subplots
        suptitle(' ') ;

        print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_atten_profile.png'] ;
        eval(print_string) ;
        pause(pause_before) ;

    % MR profile
     if (mr_empty == 0)
        if (ref_point_num == 1)
            mr_prof_mat = zeros(num_ref_points, length(range_mm)) ;
        end
        mr_prof = interp3(mr_3d, col_temp, row_temp, pln_temp) ;  % MR data
        mr_kkk = find(isnan(mr_prof)) ;
        mr_prof(mr_kkk) = 0 ;
        mr_prof_mat(ref_point_num, :) = mr_prof ;
        
        % Measurement of MR profile for inner bone to brain edge, CSF

        % Create an nx1 matrix to store values of csf_mm
        % where n = number of landmark points being examined. Initializing
        % the vector prior to calling find_edges_prof_mr to ensure 0's in case
        % find_edges_prof_mr fails.
        if (ref_point_num == 1)
           find_prof_results_mr = zeros(num_ref_points, 1) ;
        end
       
        % Find values in profile plot where 2nd derivitave crosses zero
        b_peak_i = prof_edges_new_index(bone_in1_pt); % out2?
        [mr_prof_edges_2nd, csf_edge_i, csf_mm] = ...
            find_edges_prof_mr(mr_prof, range_mm, bi_edge_i, b_peak_i, mr_3d);
        
        % store values in a matrix indexed by num_ref_points. This allows the
        % values to be written to the csv file at the end of normalcy. 
        find_prof_results_mr(ref_point_num, 1) = csf_mm ; 
        disp(['CSF (mm) = ' num2str(csf_mm)]) ;
        
        % Plot profiles of MR with CT data        
        figure(8) ;
        hold off;
        subplot(1,1,1), [AX,H1,H2] = plotyy(range_mm, prof, range_mm, mr_prof) ;
        hold on;
        axis_y_min_mm =  -1025 ;
        axis_y_max_mm = (floor(max(prof) / 100) + 1) * 100 ;
        axis(AX(1),[axis_x_min_mm axis_x_max_mm axis_y_min_mm axis_y_max_mm]) ;
        mr_axis_y_min_mm =  0 ;
        mr_axis_y_max_mm = (floor(max(mr_prof) / 100) + 1) * 100 ;
        axis(AX(2),[axis_x_min_mm axis_x_max_mm mr_axis_y_min_mm mr_axis_y_max_mm]) ;
        xlabel('Distance from Reference Point Along Surface Normal (mm)') ;
        ylabel(AX(1), 'CT Image Intensity (HU)') ;
        ylabel(AX(2), 'MR Data (Grayscale)');
      
        % Display HU levels on plots
        hold on ;
        plot([axis_x_min_mm axis_x_max_mm],[skin_thr skin_thr],'--','Color',[1 0 1]) ;
        hold on ; 
        plot([axis_x_min_mm axis_x_max_mm],[fat_thr fat_thr],'--yellow') ;
        hold on ;
        plot([axis_x_min_mm axis_x_max_mm],[mus_thr mus_thr],'--red') ;
        hold on ;
        plot([axis_x_min_mm axis_x_max_mm],[trab_thr trab_thr],'--green') ;
        hold on ;
        plot([axis_x_min_mm axis_x_max_mm],[cort_thr cort_thr],'--blue') ;             
            
        % plot profile peaks based on zero crossing of gradient                      
        hold on ;
        
        % calculate x coordinate at which to plot HU at a zero crossing (peak)
         ccc1= -10.01 + (0.01 * prof_edges_new_index) ;               
         plot (ccc1, prof_edges_new, '+yellow') ; 

        %plot first point in blue
        hold on ;
        plot (ccc1(first_pt), prof_edges_new(first_pt), '+blue') ;

        %plot second skin point in blue
        hold on ;
        plot (ccc1(skin_out2_pt), prof_edges_new(skin_out2_pt), '+blue') ;        
        
        %plot first outer bone point in cyan
        hold on ;
        plot (ccc1(bone_out1_pt), prof_edges_new(bone_out1_pt), '+red') ;

        %plot second outer bone point in cyan
        hold on ;
        plot (ccc1(bone_out2_pt), prof_edges_new(bone_out2_pt), '+red') ;        
        
        %plot first inner bone point in green
        hold on ;
        plot (ccc1(bone_in1_pt), prof_edges_new(bone_in1_pt), '+green') ; 
        
        %plot second inner bone point in green
        hold on ;
        plot (ccc1(bone_in2_pt), prof_edges_new(bone_in2_pt), '+green') ;        

        %Plot a vertical line denoting the skin edge
        hold on ;
        plot ([range_mm(sk_edge_i) range_mm(sk_edge_i)], ...
            [axis_y_min_mm axis_y_max_mm], 'blue') ;

        %Plot a vertical line denoting the outer bone edge
        hold on ;
        plot ([range_mm(bo_edge_i) range_mm(bo_edge_i)], ...
            [axis_y_min_mm axis_y_max_mm], 'red') ;

        %Plot a vertical line denoting the inner bone edge
        hold on ;
        plot ([range_mm(bi_edge_i) range_mm(bi_edge_i)], ...
            [axis_y_min_mm axis_y_max_mm], 'green') ;       
              
        %Plot a vertical line denoting the brain edge
        hold on ;
        plot ([range_mm(csf_edge_i) range_mm(csf_edge_i)], ...
            [axis_y_min_mm axis_y_max_mm], 'black') ;            
        
      title_string = ['Subject ' num2str(gray_file_name_root(1:6)) ',' ...
           '  Reference Point xyz = (' ...
           num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
           num2str(ref_point_xyz(3)) ')'] ;
       title_string1 = ['CSF Thickness (mm)= ' num2str(csf_mm)] ;
       title({title_string, title_string1}) ;

        % This plots the Main Title at the very top. Comes after all subplots
        print_string = ['print -dpng -r600 ' outputs_directory '\' ...
        out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_cr_mr_profile.png'] ;
        eval(print_string) ;
        pause(pause_before) ;
        
        %Add figure with bone inner edge and csf edge on plots
        bi_edge_rcp(1) = ref_point_rcp(1) + (range_mm(bi_edge_i)/0.5 * direction_vec_rcp(1)) ; %voxel size
        bi_edge_rcp(2) = ref_point_rcp(2) + (range_mm(bi_edge_i)/0.5 * direction_vec_rcp(2)) ;
        bi_edge_rcp(3)  = ref_point_rcp(3) + (range_mm(bi_edge_i)/0.5 * direction_vec_rcp(3)) ;
        csf_edge_rcp(1) = ref_point_rcp(1) + (range_mm(csf_edge_i)/0.5 * direction_vec_rcp(1)) ;
        csf_edge_rcp(2) = ref_point_rcp(2) + (range_mm(csf_edge_i)/0.5 * direction_vec_rcp(2)) ;
        csf_edge_rcp(3)  = ref_point_rcp(3) + (range_mm(csf_edge_i)/0.5 * direction_vec_rcp(3)) ;
        
        disp(['Displaying the transverse plane at the reference point ' ...
                'and the normal vector ...']) ;
        disp_transverse_normal_mr_csf(mr_3d, g_3d, b_3d, p_3d, ref_point_rcp, ...
                center_est_rcp, bi_edge_rcp, csf_edge_rcp, pause_before) ;
        figure(4) ;
        title_string = ['Reference Point xyz = (' ...
            num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ...
            ',' num2str(ref_point_xyz(3)) ')'] ;
        suptitle(title_string) ;
        print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_transverse.png'] ;
        eval(print_string) ;
        disp(['Displaying the sagittal plane at the reference point ' ...
            'and the normal vector ...']) ;
        disp_sagittal_normal_mr_csf(mr_3d, g_3d, b_3d, p_3d, ref_point_rcp, ...
            center_est_rcp, bi_edge_rcp, csf_edge_rcp, pause_before) ;
        figure(5) ;
        title_string = ['Reference Point xyz = (' ...
            num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
            num2str(ref_point_xyz(3)) ')'] ;
        suptitle(title_string) ;
        print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_sagittal.png'] ;
        eval(print_string) ;
        disp(['Displaying the coronal plane at the reference point ' ...
            'and the normal vector ...']) ;
        disp_coronal_normal_mr_csf(mr_3d, g_3d, b_3d, p_3d, ref_point_rcp, ...
            center_est_rcp, bi_edge_rcp, csf_edge_rcp, pause_before) ;
        figure(6) ;
        title_string = ['Reference Point xyz = (' ...
            num2str(ref_point_xyz(1)) ',' num2str(ref_point_xyz(2)) ',' ...
            num2str(ref_point_xyz(3)) ')'] ;
        suptitle(title_string) ;
        print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_coronal.png'] ;
        eval(print_string) ;
       
     end
       
        disp('*********************************************************') ;
        disp(['Name of segmentation image:         ' ...
            full_gray_file_img_name '_bin']) ;
        disp(['Name of grayscale image:            ' ...
            full_gray_file_img_name]) ;
        disp(['Radius of patch (mm) = ' num2str(r_patch_mm)]) ;
        disp(['Reference point (x,y,z) = (' ...
            num2str(ref_point_xyz(1)) ', ' ...
            num2str(ref_point_xyz(2)) ', ' ...
            num2str(ref_point_xyz(3)) ')']) ;
        disp(['Center of best-fitting sphere (x,y,z) = (' ...
            num2str(center_est_xyz(1)) ', ' ...
            num2str(center_est_xyz(2)) ', ' ...
            num2str(center_est_xyz(3)) ')']) ;
        disp(['Unit direction vector into head (x, y, z) = (' ...
            num2str(direction_vec_xyz(1)) ', ' ...
            num2str(direction_vec_xyz(2)) ', ' ...
            num2str(direction_vec_xyz(3)) ')']) ;
        disp(['Sum of squared errors when fitting sphere = ' ...
            num2str(sse_search)]) ;
        if (exitflag == 1)
            disp(['Optimization routine reported satisfactory ' ...
                  'convergence.']) ;
        else
            disp(['WARNING: OPTIMIZATION ROUTINE FAILED. EXITFLAG = ' ...
                num2str(exitflag) '.']) ;
        end
        disp('*********************************************************') ;

    end

    % Calculate the braincase volume and skull measurements
    % This next if statement allows the braincase segmentation to continue
    % and uses a default HU for the inside bone HU value. This condition
    % may occur if find_edges_prof can't find edges.
    flag = and((steps_error_vec(14) == 0), (steps_error_vec(15) == 0)) ;
    if (flag == 1)
        bi_edge_HU_values(ref_point_num) = 400;
    end
    [bii,bij,bi_values_found] = find(bi_edge_HU_values) ;
    
    % This calculates the average HU value for all interior bone edges
    % found and then takes the half-way point to use for braincase
    % segmentation.
    ave_bi_edge_hu= (mean(bi_values_found))/2 ;
    kern = ones(3,3,3) ;
    [ bc_vol, sk_l, sk_w, sk_h, sk_geo_mean, sk_perim ] = Calc_braincase...
        ( g_3d, num_row, num_col, num_pln, kern, s_3d, local_origin(3), ave_bi_edge_hu,...
          outputs_directory, out_file_name_root );
    disp(['skull_thr for braincase segmentation = ' num2str(ave_bi_edge_hu)]) ;
    disp(['Volume of braincase = ' num2str(bc_vol)]) ;
    disp(['skull length = ' num2str(sk_l)]) ;
    disp(['skull width = ' num2str(sk_w)]) ;
    disp(['skull height = ' num2str(sk_h)]) ;
    disp(['skull geometric mean = ' num2str(sk_geo_mean)]) ;
    figure(9);
    print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_bc_sagittal.png'] ;
    eval(print_string) ;
    figure(10);
    print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_bc_coronal.png'] ;
    eval(print_string) ;
    figure(11);
    print_string = ['print -dpng -r600 ' outputs_directory '\' ...
            out_file_name_root '_' num2str(ref_point_num, '%3.3d') ...
            '_bc_pfill.png'] ;
    eval(print_string) ;
   
    % Create a ".mat" file that contains all of the relevant variables and
    % matrices.

    if (mr_empty == 1)
        full_mat_file_name = [outputs_directory '\' out_file_name_root ...
            '.mat'] ;
        save(full_mat_file_name, ...
            'center_est_mat_xyz', ...
            'direction_vec_mat_xyz', ...
            'drcp_mm', ...
            'exitflag', ...
            'full_gray_file_hdr_name' , ...
            'full_gray_file_img_name', ...
            'full_ref_local_file_name', ...
            'full_ref_landmark_file_name', ...
            'gray_file_name_root', ...
            'inputs_directory', ...
            'out_file_name_root', ...
            'outputs_directory', ...
            'prof_mat', ...
            'r_patch_mm', ...
            'radius_est_vec_mm', ...
            'range_mm', ...
            'ref_point_mat_xyz', ...
            'par','pal','nas', ...
            'sse_search') ;
    else 
        full_mat_file_name = [outputs_directory '\' out_file_name_root ...
            '.mat'] ;
        save(full_mat_file_name, ...
            'center_est_mat_xyz', ...
            'direction_vec_mat_xyz', ...
            'drcp_mm', ...
            'exitflag', ...
            'full_gray_file_hdr_name' , ...
            'full_gray_file_img_name', ...
            'full_ref_local_file_name', ...
            'full_ref_landmark_file_name', ...
            'gray_file_name_root', ...
            'inputs_directory', ...
            'out_file_name_root', ...
            'outputs_directory', ...
            'full_mr_file_hdr_name', ...
            'full_mr_file_img_name', ...
            'full_mr_ref_landmark_file_name', ...
            'prof_mat', ...
            'mr_prof_mat', ...
            'r_patch_mm', ...
            'radius_est_vec_mm', ...
            'range_mm', ...
            'ref_point_mat_xyz', ...
            'par','pal','nas', ...
            'sse_search') ;
    end
    
    % Create a comma-separated-values (CSV) list which contains the
    % reference points, the estimated centers, the estimated radius, and
    % the estimated direction vectors. Added in the soft tissue thickness
    % and bone tissue thickness measures along with the average HU of the
    % bone section.

    full_csv_file_name = ...
        [outputs_directory '\' out_file_name_root '.csv'] ;
    diary(full_csv_file_name) ;

    if (mr_empty == 1)
        disp_string = ['Label, Reference x, Reference y, Reference z, ' ...
         'Estimated Center x, Estimated Center y, Estimated Center z, ' ...
         'Estimated Radius r, Estimated Direction Vector x, ' ...
         'Estimated Direction Vector y, Estimated Direction Vector z, ' ...
         'ST_mm, BT_mm, Ave_HU_Bone, BC_Vol, Skull_threshold, Skull_L, Skull_W, Skull_H, ' ...
         'Skull_perim_mm, Skull_Geo_Mean'] ;
        disp(disp_string) ;

        for ref_point_num = 1:num_ref_points

            disp_string = [ ...
            gray_file_name_root(1:9) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(radius_est_vec_mm(ref_point_num)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(find_prof_results(ref_point_num, 1)) ',' ...
            num2str(find_prof_results(ref_point_num, 2)) ',' ...
            num2str(find_prof_results(ref_point_num, 3)) ',' ...
            num2str(bc_vol) ',' ...
            num2str(ave_bi_edge_hu) ',' ...
            num2str(sk_l) ',' ...
            num2str(sk_w) ',' ...
            num2str(sk_h) ',' ...
            num2str(sk_perim) ',' ...
            num2str(sk_geo_mean)] ;
                        
            disp(disp_string) ;
        end
    else
        disp_string = ['Label, Reference x, Reference y, Reference z, ' ...
         'Estimated Center x, Estimated Center y, Estimated Center z, ' ...
         'Estimated Radius r, Estimated Direction Vector x, ' ...
         'Estimated Direction Vector y, Estimated Direction Vector z, ' ...
         'ST_mm, BT_mm, Ave_HU_Bone, CSF_mm, BC_Vol, Skull_threshold, Skull_L, Skull_W, Skull_H, ' ...
         'Skull_perim_mm, Skull_Geo_Mean, Num_Reg_Pts, SSE_Reg_vox'] ;
        disp(disp_string) ;

        for ref_point_num = 1:num_ref_points

            disp_string = [ ...
            gray_file_name_root(1:9) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(ref_point_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(center_est_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(radius_est_vec_mm(ref_point_num)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 1)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 2)) ',' ...
            num2str(direction_vec_mat_xyz(ref_point_num, 3)) ',' ...
            num2str(find_prof_results(ref_point_num, 1)) ',' ...
            num2str(find_prof_results(ref_point_num, 2)) ',' ...
            num2str(find_prof_results(ref_point_num, 3)) ',' ...
            num2str(find_prof_results_mr(ref_point_num, 1)) ',' ...
            num2str(bc_vol) ',' ...
            num2str(ave_bi_edge_hu) ',' ...
            num2str(sk_l) ',' ...
            num2str(sk_w) ',' ...
            num2str(sk_h) ',' ...
            num2str(sk_perim) ',' ...
            num2str(sk_geo_mean) ',' ...
            num2str(num_reg_pts) ',' ...
            num2str(sse_reg)] ;
                        
            disp(disp_string) ;
        end
    end

    diary off ;
    
end