%**************************************************************************
%
%  PROGRAM TITLE      find_edges_prof.m
%
%  WRITTEN BY         Kirk E. Smith
%  DATE WRITTEN       November 3, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY      
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [ prof_edges, prof_edges_new, prof_edges_new_index, first_pt...
%           skin_out2_pt, bone_out1_pt, bone_out2_pt, bone_in1_pt, bone_in2_pt ...
%           sk_edge_i, bo_edge_i, bi_edge_i, Ave_HU_Bone, ST_mm, BT_mm, steps_error_vec] ...
%           = find_edges_prof( prof, range_mm, tissue_type, tissue_thr ) ;
%
%     where
%       prof_edges          profile edges of surface normal
%       prof_edges_new      restricted length profile edges
%       prof_edges_new_index index of new edges
%       first_pt            first point on or before edge defined by binary image        
%       skin_out2_pt        point of outside of skin
%       bone_out1_pt        first point of bone outside
%       bone_out2_pt        second point of bone outside
%       bone_in1_pt         first point of bone inside
%       bone_in2_pt         second point of bone inside
%       sk_edge_i           skin edge index of profile
%       bo_edge_i           bone outside edge index of profile
%       bi_edge_i           bone inside edge index of profile
%       Ave_HU_Bone         average HU of bone
%       ST_mm               skin thickness, mm
%       BT_mm               bone thickness, mm
%       steps_error_vec     vector of errors in processing steps
%       prof                profile data of surface normal
%       range_mm            profile range, mm
%       tissue_type         tissue type vector
%       tissue_thr          tissue threshold vector
%
%  PROGRAM DESCRIPTION
%       This function examines the profile plot along the surface normal, 
%   finds egdges, and measures the thicknesses of bone and skin.
%       Changed trab_thr from 200 to 225 to account for bright skin value due to
%   kernel or other effects. Value changed in vector passed in from normalcy.   
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

function [ prof_edges, prof_edges_new, prof_edges_new_index, first_pt...
    skin_out2_pt, bone_out1_pt, bone_out2_pt, bone_in1_pt, bone_in2_pt ...
    sk_edge_i, bo_edge_i, bi_edge_i, Ave_HU_Bone, ST_mm, BT_mm, steps_error_vec] ...
    = find_edges_prof( prof, range_mm, tissue_type, tissue_thr )

prof_grad = gradient(prof, 0.01) ;
temp1 = length(prof_grad) ;
prof_grad_cross_rise = zeros(size(prof_grad)) ;
prof_grad_cross_fall = zeros(size(prof_grad)) ;
junk_type = tissue_type(1) ;
skin_type = tissue_type(2) ;
fat_type = tissue_type(3) ;
mus_type = tissue_type(4) ;
trab_type = tissue_type(5) ;
cort_type = tissue_type(6) ;
skin_thr = tissue_thr(1) ;
fat_thr = tissue_thr(2) ;
mus_thr = tissue_thr(3) ;
trab_thr = tissue_thr(4) ;
cort_thr = tissue_thr(5) ;
for index=1:temp1-1 ;
     if (prof_grad(index) >= 0);
        if (prof_grad(index+1) < 0) ;
            prof_grad_cross_rise(index) = 1 ;
        end ;
     end ;
    if (prof_grad(index) <= 0);
        if (prof_grad(index+1) > 0) ;
            prof_grad_cross_fall(index) = 1 ;
        end ;
   end ;        

end
prof_cross = prof_grad_cross_rise + prof_grad_cross_fall ;
prof_edges = prof .* prof_cross ;
% assign tissue type labels based on HU at crossings
% edit to represent sof tissue and bone only
[tyr, typc, prof_type] = find(prof_edges) ;
edges_temp = prof_type ;
[rj, cj, vj] = find(edges_temp < skin_thr) ;
prof_type(cj) = junk_type ;
[rs, cs, vs] = find(and(edges_temp >= skin_thr, edges_temp < trab_thr)) ;
prof_type(cs) = skin_type ;

[rc, cc, vc] = find(edges_temp >= trab_thr) ;
prof_type(cc) = cort_type ;

%filter edges based on classification
prof_edges_temp = zeros(size(prof_type)) ;
i_tot = length(prof_type) ;
i_end = i_tot - 1 ;
for i=1:i_end ;
    if ((prof_type(i+1) - prof_type(i)) ~= 0) ;
        prof_edges_temp(i+1) = 1 ;
    end
    if ((prof_type(i) - prof_type(i+1)) ~= 0) ;
        prof_edges_temp(i) = 1 ;
    end
    
end

%Initialize an error vector to track successful completion of steps. If a
%step is successful, then a flag will be set. The number of steps is equal
%to the number of outputs returned to normalcy. The vector can be indexed
%by the step number. 1=prof_edges_new, 2=prof_edges_new_index, 3=first_pt,...
% 4=bone_out2_pt, 5=bone_in1_pt, 6=bone_in2_pt, 7=bone_out1_pt, ...
% 8=skin_out2_pt, 9=hu_skin_edge, 10=hu_bone_out_edge, 11=hu_bone_in_edge,...
% 12=sk_edge_i, 13=bo_edge_i, 14 = bi_edge_i, 15=Ave_HU_Bone, ...
% 16=ST_mm, 17=BT_mm
steps_error_vec = zeros(1,17) ;

% Restricting range over which to search to first 20 edges to try to
% prevent cases where edge extends into another bright area. A better idea
% will be to limit the length of the profile line to begin with.
% Shortened profile length to 45 mm in normalcy
% [per, pec, pev] = find(prof_edges_temp(1:30)) ;
[per, pec, pev] = find(prof_edges_temp) ;
prof_edges_new = prof(typc(pec)) ;
%error checking
if isempty(prof_edges_new) ;
    steps_error_vec(1) = 0 ;
    disp('edge crossings not found') ;
%    return ;
else steps_error_vec(1) = 1;
end

prof_edges_new_index = typc(pec) ;
%error checking
if isempty(prof_edges_new_index) ;
    steps_error_vec(2) = 0 ;
    disp('edge crossings index not found') ;
%    return ;
else steps_error_vec(2) = 1;
end

% find first point on or before edge defined by binary image.
first_pt = max(find(prof_edges_new_index <= 1000)) ;
%error checking
if isempty(first_pt) ;
    steps_error_vec(3) = 0 ;
    disp('first point representing skin not found') ;
%    return ;
else steps_error_vec(3) = 1;
end

% Find outer and inner bone edges. This represents total bone and does not 
% look at sinus or cancellous bone between outer and inner edges. This
% finds the index into the vector that contains the edges. A problem occurs
% when skin is above the trab_thr range as bone_out2_pt then gets set; or
% if thin bone causes bone value to fall below trab_thr; Right now this is
% being trapped after return to normalcy which then breaks out of that
% landmark for loop, prints an error message, and continues on.

bone_out2_pt = min(find(prof_edges_new >= trab_thr)) ;
%error checking
if isempty(bone_out2_pt) ;
    steps_error_vec(4) = 0 ;
    disp('bone_out2_pt not found') ;
%    return ;
else steps_error_vec(4) = 1;
end

size_prof_edges_new = size(prof_edges_new);
adj_prof_edges_new = size_prof_edges_new(2);
bone_in1_pt = max(find(prof_edges_new >= trab_thr)) ;
% the following statement was added to check for the case where the last
% edge found is a bright area that is most likely due to the profile
% extending into another bone area. When that occurs, the bone_in1_pt gets
% set wrong which then causes bone_in2_pt to point outside the
% prof_edges_new vector. There is probably a better way to do this. We
% could find the max value but require it to have a brain value after it.
if (bone_in1_pt == adj_prof_edges_new)
    bone_in1_pt = max(find(prof_edges_new(1:adj_prof_edges_new-1) >= trab_thr)) ;
end
%error checking
if isempty(bone_in1_pt) ;
    steps_error_vec(5) = 0 ;
    disp('bone_in1_pt not found') ;
%    return ;
else steps_error_vec(5) = 1;
end

bone_in2_pt = bone_in1_pt +1 ;
%error checking
if isempty(bone_in2_pt) ;
    steps_error_vec(6) = 0 ;
    disp('bone_in2_pt not found') ;
%    return ;
else steps_error_vec(6) = 1;
end

bone_out1_pt = max(find(and((prof_edges_new(first_pt:bone_out2_pt) >= skin_thr),...
    (prof_edges_new(first_pt:bone_out2_pt) < trab_thr))));
%error checking
if isempty(bone_out1_pt) ;
    steps_error_vec(7) = 0 ;
    disp('bone_out1_pt not found') ;
%    return ;
else steps_error_vec(7) = 1;
end

skin_out2_pt = min(find(and((prof_edges_new(first_pt:bone_out2_pt) >= skin_thr),...
    (prof_edges_new(first_pt:bone_out2_pt) < trab_thr))));
%error checking
if isempty(skin_out2_pt) ;
    steps_error_vec(8) = 0 ;
    disp('skin_out2_pt not found') ;
%    return ;
else steps_error_vec(8) = 1;
end

% adjust index so that it points back to prof_edges_new
skin_out2_pt = first_pt -1 + skin_out2_pt ;
bone_out1_pt = first_pt -1 + bone_out1_pt ;

% find average HU for a given edge (soft tissue, bone)
hu_skin_edge=(prof_edges_new(first_pt) + prof_edges_new(skin_out2_pt))/2;
%error checking
if isempty(hu_skin_edge) ;
    steps_error_vec(9) = 0 ;
    disp('hu_skin_edge not found') ;
%    return ;
else steps_error_vec(9) = 1;
end

hu_bone_out_edge=(prof_edges_new(bone_out1_pt)+prof_edges_new(bone_out2_pt))/2;
%error checking
if isempty(hu_bone_out_edge) ;
    steps_error_vec(10) = 0 ;
    disp('hu_bone_out_edge not found') ;
%    return ;
else steps_error_vec(10) = 1;
end

hu_bone_in_edge=(prof_edges_new(bone_in1_pt)+prof_edges_new(bone_in2_pt))/2;
%error checking
if isempty(hu_bone_in_edge) ;
    steps_error_vec(11) = 0 ;
    disp('hu_bone_in_edge not found') ;
%    return ;
else steps_error_vec(11) = 1;
end

% make edge indexes point to original prof vector
sk_pt1_i = prof_edges_new_index(first_pt) ;
sk_pt2_i = prof_edges_new_index(skin_out2_pt) ;
bo_pt1_i = prof_edges_new_index(bone_out1_pt) ;
bo_pt2_i = prof_edges_new_index(bone_out2_pt) ;
bi_pt1_i = prof_edges_new_index(bone_in1_pt) ;
bi_pt2_i = prof_edges_new_index(bone_in2_pt) ;
% if any point does not exist then set all pts to 0 in order to let the
% program run through remaining landmarks
% pt0_vec = [first_pt skin_out2_pt bone_out1_pt bone_out2_pt ...
%         bone_in1_pt bone_in2_pt] ;
% if (pt0_vec(1:end) == 0) ;
%     sk_pt1_i = 0 ;
%     sk_pt2_i = 0 ;
%     bo_pt1_i = 0 ;
%     bo_pt2_i = 0 ;
%     bi_pt1_i = 0 ;
%     bi_pt2_i = 0 ;
% if the second skin point is missing, then interpolate the HU and
% search over a range 5 mm beyond the first skin point location. Skin
% interpolate at HU = 30, bone inner pt 1 = -30 (fat) Now searching over a
% range from sk_pt1 to bone outer 2. This assumes that since the bone outer
% pt 1 was not found that there is no cortical dip between outer 1 and
% outer 2
if (isempty(sk_pt2_i)) ;
    hu_skin_point = 30 ;
    hu_bone_out_point = -30 ;
%    sk_pt2_interp_i = sk_pt1_i + 500 ;
%    [skr,sk_pt2_i]=min(abs(prof(sk_pt1_i:sk_pt2_interp_i)- hu_skin_point));
    [skr,sk_pt2_i]=min(abs(prof(sk_pt1_i:bo_pt2_i)- hu_skin_point));
    sk_pt2_i = sk_pt1_i + sk_pt2_i -1 ;
    bo_pt1_i = sk_pt2_i ;
    hu_skin_edge=(prof_edges_new(first_pt) + hu_skin_point)/2;
    hu_bone_out_edge = (hu_bone_out_point+prof_edges_new(bone_out2_pt))/2;
    disp('Skin edge2 not found: Interpolating at 30 HU up to bone edge2') ;
    disp('Bone out edge not found: Interpolating at -30 HU up to bone edge2') ;
    %error checking
    if isempty(hu_skin_edge) ;
    steps_error_vec(8) = 0 ;
    disp('skin edge interpolation unsuccessful') ;
%    return ;
    else steps_error_vec(8) = 1;
         steps_error_vec(9) = 1;
    end
    if isempty(hu_bone_out_edge) ;
    steps_error_vec(10) = 0 ;
    disp('bone outer edge interpolation unsuccessful') ;
%    return ;
    else steps_error_vec(7) = 1;
         steps_error_vec(10) = 1;
    end
    if (isempty(bo_pt1_i)) ;
         bo_pt1_i = sk_pt2_i ;
         hu_bone_out_edge = (hu_bone_out_point+prof_edges_new(bone_out2_pt))/2;
         disp('Bone out edge not found: Interpolating at -30 HU up to bone edge2') ;
    %error checking
        if isempty(hu_bone_out_edge) ;
        steps_error_vec(10) = 0 ;
        disp('bone outer edge interpolation unsuccessful') ;
%    return ;
        else steps_error_vec(7) = 1;
         steps_error_vec(10) = 1;
        end
    end     
end
% find index where average HU for an edge is closest to that location
[skr,sk_edge_i]=min(abs(prof(sk_pt1_i:sk_pt2_i)- hu_skin_edge));
[bor,bo_edge_i]=min(abs(prof(bo_pt1_i:bo_pt2_i)- hu_bone_out_edge));
[bor,bi_edge_i]=min(abs(prof(bi_pt1_i:bi_pt2_i)- hu_bone_in_edge));

% adjust index so that it points back to prof
sk_edge_i = sk_pt1_i + sk_edge_i -1 ;
%error checking
if isempty(sk_edge_i) ;
    steps_error_vec(12) = 0 ;
    disp('sk_edge_i not found') ;
%    return ;
else steps_error_vec(12) = 1;
end

bo_edge_i = bo_pt1_i + bo_edge_i -1 ;
%error checking
if isempty(bo_edge_i) ;
    steps_error_vec(13) = 0 ;
    disp('bo_edge_i not found') ;
%    return ;
else steps_error_vec(13) = 1;
end

bi_edge_i = bi_pt1_i + bi_edge_i - 1 ;
%error checking
if isempty(bi_edge_i) ;
    steps_error_vec(14) = 0 ;
    disp('bi_edge_i not found') ;
%    return ;
else steps_error_vec(14) = 1;
end

% calculate average HU value over region
tot_bone_pts = length(prof(bo_edge_i:bi_edge_i)) ;
Ave_HU_Bone = (sum(prof(bo_edge_i:bi_edge_i))/tot_bone_pts) ;
%error checking
if or((isempty(Ave_HU_Bone)), (isnan(Ave_HU_Bone))) ;
    steps_error_vec(15) = 0 ;
    disp('Ave_HU_Bone not found') ;
%    return ;
else steps_error_vec(15) = 1;
end

% calculate soft tissue thickness and bone tissue thickness in mm
ST_mm = range_mm(bo_edge_i) - range_mm(sk_edge_i) ;
%error checking
if isempty(ST_mm) ;
    steps_error_vec(16) = 0 ;
    disp('ST_mm not found') ;
%    return ;
else steps_error_vec(16) = 1;
end

BT_mm = range_mm(bi_edge_i) - range_mm(bo_edge_i) ;
%error checking
if isempty(BT_mm) ;
    steps_error_vec(17) = 0 ;
    disp('BT_mm not found') ;
%    return ;
else steps_error_vec(17) = 1;
end

end

