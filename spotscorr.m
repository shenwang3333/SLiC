function [out] = spotscorr (img1, img2, res_int1, cvt_flag, disp_flag)

% spotscorr: spots correlation designed for visualization of peptide/protein-liposome binding assay,
%			 generally for SLiC assays, but may be also adaptive for similar schemes. The input requires
%			 a calculated output from the function spotdet or spotmulsz, which is the detection result 
%			 of img1; img2 is the same vision compared to img1, the difference between them is the 
%			 fluorescence channel, e.g., img1 is green, img2 is red/far-red, and they are the same 
%			 vision. spotscorr will display the two images with detections calculated in img1.

% INPUT parameters:
%	img1: 2d single channel image with a certain type, the detected image
%	img2: 2d single channel image with a certain type, with a same vision to img1
%	res_int1: n-by-4 array, detection results of img1 use spotdet or spotmulsz, the columns are row-
%			  coordinates, col-coordinates, integrated spot intensities, and object scale flags
%	cvt_flag: convertion flag. Before using this flag, one should measure the liposome size distribution
%			  by dynamic light scattering, then imaging the same liposome to aquire itensity distribution.
%			  First manually do natural logarithm processing for the two datasets, then fitting the 
%			  distribution with Gaussian or Gaussian Amplitude function. If cvt_flag is set to a non-zero 
%			  value, spotscorr will ask user to type in two values, one is the expectation value of the 
%			  size distribution (natural logarithm processed) measured by DLS, the other is the expectation 
%			  value of the intensity profile distribution (natural logarithm processed) measured by 
%			  fluorescence schemes. spotscorr will convert the intensity profiles of liposomes to real sizes 
%			  that are measured by dynamic light scattering
%	disp_flag: display flag. If it is set to a non-zero value, spotscorr will plot linear and log style 
%			   I(liposome)/Size(liposome)-I(peptide) sensing graph

% OUTPUT parameters:
%	out: n-by-2 array contains the integrated intensities/liposome sizes of img2 that echos img1's 
%		 detection (first column) and integrated intensities of the spots in img1 (second column)

% Written by Shen Wang, Sep. 19th, 2018, in HUST



img2 = double(img2);

[nspot, ~] = size(res_int1);
flags_col = res_int1(:, 4);

%% find unique object scale in res_int1
unique_flags = unique(flags_col);
max_flags = max(unique_flags);

%% create cell array to accomodate masks, object scale value as the index of cell array
CC_masks = cell(max_flags, 1);
cc_mask_sur = zeros(max_flags, 1);
len_unique_flags = length(unique_flags);

for i = 1:len_unique_flags

	mask = zeros(unique_flags(i));
	r = floor(unique_flags(i)/2);
	mask(r+1, r+1) = 1;
	mask_bw = im2bw(mask);
	SE = strel('disk', r, 0);
	mask_bw = imdilate(mask_bw, SE);
	mask = double(mask_bw);
	CC_masks{unique_flags(i)} = mask;
	cc_mask_sur(unique_flags(i)) = sum(sum(mask));

end


%% calculate img2 integrated spot intensities which are correlated to the spots detected in img1
res_int2 = zeros(nspot, 1);
os_sur = zeros(nspot, 1);

for j = 1:nspot

	spot_window = img2(res_int1(j, 1)-floor(flags_col(j)/2):res_int1(j, 1)+floor(flags_col(j)/2), ...
		res_int1(j, 2)-floor(flags_col(j)/2):res_int1(j, 2)+floor(flags_col(j)/2));
	spot_window_pmsk = spot_window.*CC_masks{flags_col(j)};
	int_den = sum(sum(spot_window_pmsk));
	res_int2(j) = int_den;
	os_sur(j) = cc_mask_sur(flags_col(j))

end


%% convert intensity profiles to liposome sizes
 
if cvt_flag 

	while 1

		raw_ctr = input(['please type in the expectation value of the size distribution' ...
		' ' '(natural logarithm processed) measured by DLS : ']);
		
		if raw_ctr <= 0 || ~isnumeric(raw_ctr)
			disp('input value should be a positive number, please retype.');
			continue;
		else
			dls_model_mean = raw_ctr;
			break;
		end
	
	end
	
	while 1
		
		raw_ctr2 = input(['please type in the expectation value of the intensity profile distribution'...
		 ' ' '(natural logarithm processed) measured by fluorescence schemes : ']);
		
		if raw_ctr2 <= 0 || ~isnumeric(raw_ctr2)
			disp('input center value should be a positive number, please retype.');
			continue;
		else
			intensity_mean = raw_ctr2;
			break;
		end

	end

	offset_pow = dls_model_mean/intensity_mean;
	res_int2 = res_int2.^offset_pow;

end

corrlst = zeros(nspot, 2);
corrlst(:, 2) = res_int1(:, 3)./os_sur;
corrlst(:, 1) = res_int2;


%% display linear and log style I(liposome)/Size(liposome)-I(peptide) sensing graph
if disp_flag

	figure;
	plot(corrlst(:, 1), corrlst(:, 2), '.');
	title('Linear style');

	if cvt_flag
		xlabel('Size (nm)');
		ylabel('Binding densities (a.u.)');
	else
		xlabel('Liposome densities (a.u.)');
		ylabel('Binding densities (a.u.)');		
	end

	hold on;
	
	figure;
	loglog(corrlst(:, 1), corrlst(:, 2), '.');
	title('Log style');

	if cvt_flag
		xlabel('Size (nm)');
		ylabel('Binding densities (a.u.)');
	else
		xlabel('Liposome densities (a.u.)');
		ylabel('Binding densities (a.u.)');
	end

	hold on;

end

%% display peptide detections and liposome echos
img1 = double(img1);
dispdet(img1, res_int1);
title('Peptide detections');
res_int2_complete = res_int1;
res_int2_complete(:, 3) = res_int2;
dispdet(img2, res_int2_complete);
title('Liposome echos');

out = corrlst;




