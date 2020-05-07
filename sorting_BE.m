% Copyright (C) 2020 Andreas Bertsatos <abertsatos@biol.uoa.gr>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
%
function [varargout] = sorting_BE(descriptives, left_side, right_side)
	% function sorted = sorting("descriptives.csv", "left_side.csv", "right_side.csv")
	% function [sorted, stats] = sorting("descriptives.csv", "left_side.csv", "right_side.csv")
	% function [sorted, stats, unsorted] = sorting("descriptives.csv", "left_side.csv", "right_side.csv")
  %
  % This function applies a sorting algorithm on left and right side skeletal elements based on the
	% variables produced by the GNU Octave CSG Toolkit according to the descriptives provided in the
	% relevant file.
	%
	% The csv files related to the left and right side elements should be (N+1)x47 in size, where N is
	% the number of individual bones. The first collumn should be a numeric index followed by the same
	% order of column variables found in the descriptives csv file and the first row is the header
	% with the variable names.
	%
	% The function prints a summary of the sorted elements and it returns to the user a variable number
	% of output arguments. If one output argument is given, then the sorted elements are returned in
	% a Nx3 matrix, where N is the number of sorted individuals, the first column corresponds to left
	% side, the second column to the right side elements and the third column to the related sorting
	% scores. If the sorting was accomplished through mutual exclusion, then the score value is NaN.
	% Likewise, for single elements that belong to distinct individuals, only one column contains the
	% element's numeric index and the complementary side is set to NaN along with the score (3rd column)
	%
	% When two output arguments are given, then the function will further return a 2x7 cell array
	% with information about the sorted elements and the identified individuals as well as the
	% numbers of remaining plausible pairs and definite mismatched pairs.
	% The first row is the description header for each value as shown below:
	%
	% 		{"sample size", "sorted elements", "paired elements", "single elements",...
	%			 "individuals", "remaining plausible pairs", "definite mismatched pairs"}
	%
	% If a third output argument is provided, the function will return the remaining (if any) list of
	% plausible matches that remain unsorted. The unsorted elements are also returned in a Nx3 matrix,
	% where N is the number of unsorted plausible matches. If no such cases exist, then an empty
	% matrix is returned.
	
	% check the number of input arguments
  if nargin != 3
    printf("Invalid number of input arguments. Type 'help sorting_BE' for more info\n");
    return;
  endif
  % check all imput arguments are strings
  if !ischar(descriptives) || !ischar(left_side) || !ischar(right_side')
    printf("All input arguments should be char strings of relevant csv files. Type 'help sorting_BE' for info\n");
    return;
  endif
	% check the number of output arguments
	if nargout < 1 || nargout > 3
    printf("Invalid number of output arguments. Type 'help sorting_BE' for more info\n");
    return;
  endif
	
	% load bone descriptives from file
	BD = csvread(descriptives);
	mean = BD(2,[2:end]);
	StdDev = BD(3,[2:end]);
	lbound = BD(4,[2:end]);
	ubound = BD(5,[2:end]);
	
	% load left and right side element data
	left_sample = csvread(left_side);
	left_sample(1,:) = [];
	left_sample_list = left_sample(:,1);
	left_sample_size = length(left_sample_list);
	right_sample = csvread(right_side);
	right_sample(1,:)= [];
	right_sample_list = right_sample(:,1);
	right_sample_size = length(right_sample_list);
	sample_size = left_sample_size + right_sample_size;
	
	% create a testing matrix
	TestDATA = suffle(left_sample, right_sample);
	
	% for every individual variable find the definite mismatches and append them to the reject list
	rejected = [];
	for var=1:length(lbound)
		rej_indx = unique(find(TestDATA(:,var+2) < lbound(var)));
		rejected = [rejected; rej_indx];
		rej_indx = unique(find(TestDATA(:,var+2) > ubound(var)));
		rejected = [rejected; rej_indx];
		rejected = unique(rejected);
	endfor
	
	% remove all rejected pairs from the data set and make a list of plausible matches
	plausible = TestDATA;
	plausible(rejected,:) = [];
	definite_mismatched_pairs = length(rejected);
	
	% if all possible matches have been rejected, create the sorted list according to the input samples
	if isempty(plausible)
		sorted = [];
		unsorted_pairs = 0;
		for i=1:length(left_sample_list)
			sorted = [sorted; left_sample_list(i), NaN, NaN];
		endfor
		for i=1:length(right_sample_list)
			sorted = [sorted; NaN, right_sample_list(i), NaN];
		endfor
		% calculate the number of single elements sorted by mutual exclusion
		sorted_by_elimination = size(sorted,1);
		% and return the appropriate number of output arguments
		if nargout == 1
			varargout{1} = sorted;
		elseif nargout == 2
			varargout{1} = sorted;
			varargout{2} = {"sample size", "sorted elements", "paired elements", "single elements",...
											"individuals", "unsorted_pairs", "definite_mismatched_pairs";...
											sample_size, sorted_elements, sorted_pairs, sorted_by_elimination,...
											individuals, unsorted_pairs, definite_mismatched_pairs};
		elseif nargout == 3
			varargout{1} = sorted;
			varargout{2} = {"sample size", "sorted elements", "paired elements", "single elements",...
											"individuals", "unsorted_pairs", "definite_mismatched_pairs";...
											sample_size, sorted_elements, sorted_pairs, sorted_by_elimination,...
											individuals, unsorted_pairs, definite_mismatched_pairs};
			varargout{3} = [];
		endif
		return
	endif
	
	% compare the initial left and right side lists of elements against the remaining plausible matches
	% to find potential samples belonging to individuals, who are represented by a single side element
	index = 0;
	for i=1:length(plausible(:,1))
		left_sample_list(left_sample_list == plausible(i,1)) = [];
		right_sample_list(right_sample_list == plausible(i,2)) = [];
	endfor
	sorted = [];
	for i=1:length(left_sample_list)
		sorted = [sorted; left_sample_list(i), NaN, NaN];
	endfor
	for i=1:length(right_sample_list)
		sorted = [sorted; NaN, right_sample_list(i), NaN];
	endfor
	% calculate the number of single elements sorted by mutual exclusion
	sorted_by_elimination = size(sorted,1);
	% calculate the sum of absolute z-scores for each paired samples in the testing pool
	testing = [plausible(:,[1:2]), sum((abs((plausible(:,[3:end]) - mean) ./ StdDev)), 2)];

	% scan through the remaining testing cases and cluster the associated pairs into separate subgroups
	% for each side of bones
	cluster = cluster_pairs(testing);
	
	% compare the scores of every element with each paired association between both sides and keep
	% the matching pairs as a sorted pair when both sides exhibit the lowest score on the same matched
	% pair and the score is progressively below 30
	[sorted, cluster] = compare_scores_with_threshold(sorted, cluster);
	
	% calculate the number of pairs sorted by lowest score below 30
	sorted_by_score_threshold = size(sorted,1) - sorted_by_elimination;
	
	% check if cluster exists and also contains any remaining plausible matches
	if exist('cluster') && !isempty(cluster)
		s = 0;
		for c=1:length(cluster)
			s += size(cluster(c).left,1) + size(cluster(c).right,1);
		endfor
		if s > 0
			% scan through the remaining testing cases and cluster the associated pairs into separate
			% subgroups for each side of bones
			cluster = recluster_pairs(cluster);
	
			% in each subgroup compare the scores of every element with each paired association between
			% both sides and keep the matching pairs as a sorted pair when both sides exhibit the lowest
			% score on the same matched pair and the score is lower than the second smaller one by at
			% least 5 units
			[sorted, cluster] = compare_scores_with_difference(sorted, cluster);
	
			% calculate the number of pairs sorted by lowest score with difference above 5
			sorted_by_score_difference = size(sorted,1) - (sorted_by_elimination + sorted_by_score_threshold);
			% calculate the total number of pair sorted
			sorted_pairs = sorted_by_score_threshold + sorted_by_score_difference;
		else
			sorted_pairs = sorted_by_score_threshold;
		endif
	else
		sorted_pairs = sorted_by_score_threshold;
	endif
	
	% check if cluster exists and also contains any remaining plausible matches and concatenate
	% the remaining unsorted pairs in a single matrix
	if exist('cluster') && !isempty(cluster)
		unsorted = [];
		for c=1:length(cluster)
			unsorted = [unsorted; cluster(c).left; cluster(c).right];
		endfor
		unsorted = unique(unsorted, "rows");
	else
		unsorted = [];
	endif
	unsorted_pairs = size(unsorted,1);
		
	% calculate the total number of sorted elements and distinct individuals identified
	sorted_elements = sum(sum(isfinite(sorted(:,[1:2]))));
	individuals = size(sorted,1);
	
	% report the sorting performance statistics
	printf("%i out of %i elements have been sorted into %i individuals.\n",...
					sorted_elements, sample_size, individuals);
	
	% return the appropriate number of output arguments
	if nargout == 1
		varargout{1} = sorted;
	elseif nargout == 2
		varargout{1} = sorted;
		varargout{2} = {"sample size", "sorted elements", "paired elements", "single elements",...
											"individuals", "unsorted_pairs", "definite_mismatched_pairs";...
											sample_size, sorted_elements, sorted_pairs, sorted_by_elimination,...
											individuals, unsorted_pairs, definite_mismatched_pairs};
	elseif nargout == 3
		varargout{1} = sorted;
		varargout{2} = {"sample size", "sorted elements", "paired elements", "single elements",...
											"individuals", "unsorted_pairs", "definite_mismatched_pairs";...
											sample_size, sorted_elements, sorted_pairs, sorted_by_elimination,...
											individuals, unsorted_pairs, definite_mismatched_pairs};
		varargout{3} = unsorted;
	endif
endfunction


function TestDATA = suffle(A, B)
	% scan through both matrices and create a matrix with all possible pairs between left and right
	% side. Return a testing data set with the first column containing the index of the left side and
	% the second column containing the right side.
	index = 0;
	for i=1:length(A(:,1))
		for k=1:length(B(:,1))
			index += 1;
			% for each possible match extract the respective variable vector from each side, calculate
			% the element-wise difference of absolute values and append it in the returning testing data
			Lvec = abs(A(i,[2:end]));
			Rvec = abs(B(k,[2:end]));
			match = Lvec - Rvec;
			TestDATA(index,:) = [A(i), B(k), match];
		endfor
	endfor
endfunction

function	cluster = cluster_pairs(testing);
	% scan through the plausible paired matches and cluster the associated pairs
	% into separate subgroups for each side of bones
	savelist = testing;
	group = 0;
	while (length(testing(:,1)) > 0)
		group += 1; cluster(group).left = [];
		complete = false;
		% find a sample with minimum occurence and use it as a seed
		samples = unique(testing(:,1));
		clear nsamples;
		for s=1:length(samples)
			idx = find(testing(:,1) == samples(s));
			nsamples(s,:) = [length(idx), samples(s)];
		endfor
		nsamples = sortrows(nsamples, 1);
		idx = find(testing(:,1) == nsamples(1,2));	
		left_seed(group) = testing(idx(1),1);
		right_seed(group) = testing(idx(1),2);
		left_samples = left_seed(group);
		right_samples = [];
		while (!complete)
			% find occurences of right samples according to the left samples
			for i=1:length(left_samples)
				idx = find(testing(:,1) == left_samples(i));
				right_samples = [right_samples; testing(idx,2)];
				cluster(group).left = [cluster(group).left; testing(idx,:)];
				testing(idx,:) = [];
			endfor
			k = 0;
			% find occurences of left samples according to the right samples
			for i=1:length(right_samples)
				idx = find(testing(:,2) == right_samples(i));
				if !isempty(idx)
					left_samples = [left_samples; testing(idx,1)];
				else
					k += 1;
				endif
			endfor
			if i == k
				complete = true;
			endif
		endwhile
	endwhile
	testing = savelist;
	group = 0;
	while (length(testing(:,2)) > 0)
		group += 1; cluster(group).right = [];
		complete = false;
		right_samples = right_seed(group);  %testing(1,3);
		left_samples = [];
		while (!complete)
			% find occurences of left samples according to the right samples
			for i=1:length(right_samples)
				idx = find(testing(:,2) == right_samples(i));
				left_samples = [left_samples; testing(idx,1)];
				cluster(group).right = [cluster(group).right; testing(idx,:)];
				testing(idx,:) = [];
			endfor
			k = 0;
			% find occurences of right samples according to the left samples
			for i=1:length(left_samples)
				idx = find(testing(:,1) == left_samples(i));
				if !isempty(idx)
					right_samples = [right_samples; testing(idx,2)];
				else
					k += 1;
				endif
			endfor
			if i == k
				complete = true;
			endif
		endwhile
	endwhile
endfunction

function [sorted, cluster] = compare_scores_with_threshold(sorted, cluster);
	% in each subgroup compare the scores of every element with each paired association between
	% both sides and keep the matching pairs as a sorted pair when both sides exhibit the lowest
	% score on the same matched pair and the score is progressively below 30
	for range=20:1:30
		for c=1:length(cluster)
			% make a list of unique elements
			left_list = unique(cluster(c).left(:,1));
			for s=1:length(left_list)
				idx = find(cluster(c).left(:,1) == left_list(s));
				if length(idx) > 0
					left_side = cluster(c).left(idx,:);
					left_side = sortrows(left_side, 3);
					right_sample = left_side(1,2);
					idx = find(cluster(c).right(:,2) == right_sample);
					right_side = cluster(c).right(idx,:);
					right_side = sortrows(right_side, 3);
					score = right_side(1,3);
					% check if elements in first rows (lowest scores) match (have the same sample indices)
					% if they match append the pair in the sorted list
					if score < range && left_side(1,1) == right_side(1,1) && left_side(1,2) == right_side(1,2)
						sorted = [sorted; left_side(1,:)];
						% exclude sorted elements from the cluster
						idx = find(cluster(c).left(:,1) == left_side(1,1));
						cluster(c).left(idx,:) = [];
						idx = find(cluster(c).left(:,2) == left_side(1,2));
						cluster(c).left(idx,:) = [];
						idx = find(cluster(c).right(:,2) == right_side(1,2));
						cluster(c).right(idx,:) = [];
						idx = find(cluster(c).right(:,1) == right_side(1,1));
						cluster(c).right(idx,:) = [];
					endif
				endif
			endfor
		endfor
		% check if any subgroup contains an identical single match and in such case consider it true match
		% and append it in the sorted list and remove it from the cluster
		for c=length(cluster):-1:1
			if length(cluster(c).left(:,1)) == 1 && length(cluster(c).right(:,1)) == 1
				left_side = cluster(c).left(1,:);
				right_side = cluster(c).right(1,:);
				if left_side(1,1) == right_side(1,1) && left_side(1,2) == right_side(1,2)
					sorted = [sorted; left_side(1,:)];
					cluster(c) = [];
				endif
			else
				index = 0;
				for s=1:length(cluster(c).left(:,1))
					left_sample = cluster(c).left(s,1);
					right_sample = cluster(c).left(s,2);
					if sum(cluster(c).left(:,1) == left_sample) == 1 && sum(cluster(c).left(:,2) == right_sample) == 1
						index += 1;
						sorted = [sorted; cluster(c).left(s,:)];
						remove(index) = s;
						% remove from right cluster as well
						idx = find(cluster(c).right(:,1) == left_sample);
						cluster(c).right(idx,:) = [];
						idx = find(cluster(c).right(:,2) == right_sample);
						cluster(c).right(idx,:) = [];
					endif	
				endfor
				if exist("remove", "var")
					cluster(c).left(remove,:) = []; clear remove;
				endif
			endif
		endfor
	endfor
endfunction

function cluster2 = recluster_pairs(cluster)
	% scan through the remaining testing cases and cluster the associated pairs into separate subgroups
	% for each side of bones
	Lgroup = 0; Rgroup = 0;
	for c=1:length(cluster)
		testing = cluster(c).left;
		% scan through the remaining testing cases and cluster the associated pairs into separate subgroups
		% for each side of bones
		while (length(testing(:,1)) > 0)
			Lgroup += 1; cluster2(Lgroup).left = [];
			complete = false;
			% find a sample with minimum occurence and use it as a seed
			samples = unique(testing(:,1));
			clear nsamples;
			for s=1:length(samples)
				idx = find(testing(:,1) == samples(s));
				nsamples(s,:) = [length(idx), samples(s)];
			endfor
			nsamples = sortrows(nsamples, 1);
			idx = find(testing(:,1) == nsamples(1,2));	
			left_seed(Lgroup) = testing(idx(1),1);
			right_seed(Lgroup) = testing(idx(1),2);
			left_samples = left_seed(Lgroup);
			right_samples = [];
			while (!complete)
				% find occurences of right samples according to the left samples
				for i=1:length(left_samples)
					idx = find(testing(:,1) == left_samples(i));
					right_samples = [right_samples; testing(idx,2)];
					cluster2(Lgroup).left = [cluster2(Lgroup).left; testing(idx,:)];
					testing(idx,:) = [];
				endfor
				k = 0;
				% find occurences of left samples according to the right samples
				for i=1:length(right_samples)
					idx = find(testing(:,2) == right_samples(i));
					if !isempty(idx)
						left_samples = [left_samples; testing(idx,1)];
					else
						k += 1;
					endif
				endfor
				if i == k
					complete = true;
				endif
			endwhile
		endwhile
		testing = cluster(c).right;
		while (length(testing(:,2)) > 0)
			Rgroup += 1; cluster2(Rgroup).right = [];
			complete = false;
			right_samples = right_seed(Rgroup);
			left_samples = [];
			while (!complete)
				% find occurences of left samples according to the right samples
				for i=1:length(right_samples)
					idx = find(testing(:,2) == right_samples(i));
					left_samples = [left_samples; testing(idx,1)];
					cluster2(Rgroup).right = [cluster2(Rgroup).right; testing(idx,:)];
					testing(idx,:) = [];
				endfor
				k = 0;
				% find occurences of right samples according to the left samples
				for i=1:length(left_samples)
					idx = find(testing(:,1) == left_samples(i));
					if !isempty(idx)
						right_samples = [right_samples; testing(idx,2)];
					else
						k += 1;
					endif
				endfor
				if i == k
					complete = true;
				endif
			endwhile
		endwhile
	endfor
endfunction

function [sorted, cluster] = compare_scores_with_difference(sorted, cluster)
	% in each subgroup compare the scores of every element with each paired association between both
	% sides and keep the matching pairs as a sorted pair when both sides exhibit the lowest score
	% on the same matched pair and the score is lower than the second smaller one by at least 5 units
	for t=1:2
		for c=length(cluster):-1:1
			% make a list of unique elements
			left_list = unique(cluster(c).left(:,1));
			for s=1:length(left_list)
				idxL = find(cluster(c).left(:,1) == left_list(s));
				idxR = find(cluster(c).right(:,2) == left_list(s));
				if length(idxL) > 0 && length(idxR) > 0
					left_side = cluster(c).left(idxL,:);
					right_side = cluster(c).right(idxR,:);
					% if multiple pairs are present on right side only, use it explicitly
					if length(right_side(:,1)) > 1 && length(left_side(:,1)) == 1
						left_side = sortrows(left_side, 3);
						right_side = sortrows(right_side, 3);
						score_diff = abs(right_side(1,3) - right_side(2,3));
						if score_diff > 5 && left_side(1,1) == right_side(1,1) && left_side(1,2) == right_side(1,2)
							sorted = [sorted; left_side(1,:)];
							% exclude sorted elements from the cluster
							idx = find(cluster(c).left(:,1) == left_side(1,1));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).left(:,2) == left_side(1,2));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).right(:,2) == right_side(1,2));
							cluster(c).right(idx,:) = [];
							idx = find(cluster(c).right(:,1) == right_side(1,1));
							cluster(c).right(idx,:) = [];
						endif
					% if multiple pairs are present on left side only, use it explicitly
					elseif length(left_side(:,1)) > 1 && length(right_side(:,1)) == 1
						left_side = sortrows(left_side, 3);
						right_side = sortrows(right_side, 3);
						score_diff = abs(left_side(1,3) - left_side(2,3));
						if score_diff > 5 && left_side(1,1) == right_side(1,1) && left_side(1,2) == right_side(1,2)
							sorted = [sorted; left_side(1,:)];
							% exclude sorted elements from the cluster
							idx = find(cluster(c).left(:,1) == left_side(1,1));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).left(:,2) == left_side(1,2));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).right(:,2) == right_side(1,2));
							cluster(c).right(idx,:) = [];
							idx = find(cluster(c).right(:,1) == right_side(1,1));
							cluster(c).right(idx,:) = [];
						endif
					% if multiple pairs are present on both sides, find the minimum difference from either side
					elseif length(left_side(:,1)) > 1 && length(right_side(:,1)) > 1
						left_side = sortrows(left_side, 3);
						right_side = sortrows(right_side, 3);
						score_L = abs(left_side(1,3) - left_side(2,3));
						score_R = abs(right_side(1,3) - right_side(2,3));
						score_diff = min([score_L, score_R]);
						if score_diff > 5 && left_side(1,1) == right_side(1,1) && left_side(1,2) == right_side(1,2)
							sorted = [sorted; left_side(1,:)];
							% exclude sorted elements from the cluster
							idx = find(cluster(c).left(:,1) == left_side(1,1));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).left(:,2) == left_side(1,2));
							cluster(c).left(idx,:) = [];
							idx = find(cluster(c).right(:,2) == right_side(1,2));
							cluster(c).right(idx,:) = [];
							idx = find(cluster(c).right(:,1) == right_side(1,1));
							cluster(c).right(idx,:) = [];
						endif
					else
						sorted = [sorted; left_side(1,:)];
						% exclude sorted elements from the cluster
						idx = find(cluster(c).left(:,1) == left_side(1,1));
						cluster(c).left(idx,:) = [];
						idx = find(cluster(c).left(:,2) == left_side(1,2));
						cluster(c).left(idx,:) = [];
						idx = find(cluster(c).right(:,2) == right_side(1,2));
						cluster(c).right(idx,:) = [];
						idx = find(cluster(c).right(:,1) == right_side(1,1));
						cluster(c).right(idx,:) = [];
					endif
				endif
			endfor
		endfor
	endfor	
endfunction