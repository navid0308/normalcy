%**************************************************************************
%
%  PROGRAM TITLE      find_edges_prof_mr.m
%
%  WRITTEN BY         Gregroy G. Reiker
%  DATE WRITTEN       August 25, 2011
%  WRITTEN FOR        Pediatric head modeling project
%
%  REVISIONS BY       Gregory G. Reiker  
%  DATE MODIFIED      September 1, 2011
%  DATE MODIFIED      November 3, 2011
%
%  CALLING SYNTAX
%     Use the following syntax:
%        [mr_prof_edges_2nd, csf_edge_i, csf_mm] = ...
%           find_edges_prof_mr(mr_prof, range_mm, bi_edge_i, b_peak_i, mr_3d)  ;
%
%     where
%       mr_prof_edges_2nd   MR profile edges of second derivative
%       csf_edge_i          CSF edge index
%       csf_mm              CSF thickness, mm
%       mr_prof             MR profile of surface normal
%       range_mm            profile range, mm
%       bi_edge_i           bone inside edge index
%       b_peak_i            bonde peak index
%       mr_3d               MR volume data
%
%  PROGRAM DESCRIPTION
%   This function examines the MR profile plot along the surface normal, 
%   finds brain edge, and measures the CSF thickness.  
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
%
%
%  VERSION HISTORY
%     Version      Date                          Comment
%     -------   ---------------    ----------------------------------------
%       1.0     August 25, 2011   Initial release.
%       1.1     November 3, 2011  Added moving average of profile.
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

function [mr_prof_edges_2nd, csf_edge_i, csf_mm] = ...
            find_edges_prof_mr(mr_prof, range_mm, bi_edge_i, b_peak_i, mr_3d) 

% moving ave.
win=ones(5,1)/5;
mr_prof = conv(mr_prof,win,'same');

prof_grad = gradient(gradient(mr_prof, 0.01)) ;  % 2nd derivative
temp1 = length(prof_grad) ;
prof_grad_cross_rise = zeros(size(prof_grad)) ;
prof_grad_cross_fall = zeros(size(prof_grad)) ;


   for index=1:temp1-1-3 ;  %maybe another -3 due to filter???
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
mr_prof_edges_2nd = mr_prof .* prof_cross ;

% find brain edge and distance from inner bone edge

% possibly use to verify brain edge
%hi=max(max(max(mr_3d)));
%brain_thr=round(graythresh(mr_3d/hi)*hi); % brain threshold 

% find min MR profile value after peak of CT bone
% could use an if statement here. if b_peak_i+1000 > length(mr_prof) then
% set range to be length(mr_prof) or else set error message.
if b_peak_i+1000 > length(mr_prof)
    search_end = 5 ;
else    
search_end = 1000;
end
    for i=1:search_end
        if b_peak_i+search_end > length(mr_prof)-2 %-2 due to 5 window filter
        search_end = search_end-1;
        end
    end
[mr_bone_min, mr_bone_min_i] = min(mr_prof(b_peak_i:(b_peak_i+search_end))); %10mm/0.01 delta
mr_bone_min_i=mr_bone_min_i+(b_peak_i-1); %back to orig index ref.

% 1st crossing after inner bone edge and peak of CT bone
ind=find(prof_grad_cross_rise);
i=find((ind > bi_edge_i) & (ind >= mr_bone_min_i),1); 
csf_edge_i=ind(i);
csf_mm = range_mm(csf_edge_i)- range_mm(bi_edge_i);
end

