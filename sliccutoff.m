function [out] = sliccutoff (res_raw, lwsz_flag, nsp_flag)

% sliccutoff: automatic data cutoff designed for SLiC assay. sliccutoff have two optional flag, one is lower-size
% 			  cutoff for liposome size, a global cutoff that is independent of densities/intensities; the other
% 			  is non-specific binding cutoff, which removes echo points in an area that lower than a proper size
% 			  and lower than a proper density/intensity. One should carefully evaluate the detections and echoes,
% 			  mostly, the non-specific detections and echoes may be made by pixel noise or over-filtration, alth- 
%			  ough the non-specific echoes are sometimes avoidless.

% INPUT parameters:
%
% 	res_raw: n-by-2 array, usually the output made by spotscorr. The first column is liposome sizes or intensities
% 			 and the second column is peptide binding densities.
%	lwsz_flag: lower cutoff flag, if it is set to a non-zero value, sliccutoff will ask user to type in a value,
% 			   echoes below the value will be removed.
% 	nsp_flag: non-specific detections and echoes cutoff flag, if it is set to a non-zero value, sliccutoff will
% 			  ask user to type in two values, respectively, one is liposome size, the other is peptide binding
% 			  density, the two values make up of a rectangle region on size-density graph, any echo points in the
% 			  rectangle region will be removed.

% OUTPUT parameter:
%
% 	out: n-by-2 array, cutoff result of input made by spotscorr.

% Written by Shen Wang, Sep. 30th, 2018, in HUST


if ~lwsz_flag & ~nsp_flag

	out = res_raw;
	disp('NOTE: raw data were not rescaled/cutted. ');

end

if lwsz_flag

	while 1

		raw_lwsz = input('please type in lower cutoff of liposome size: ');

		if raw_lwsz <= 0 || ~isnumeric(raw_lwsz)
			disp('input value should be a positive number, please retype.');
			continue;
		else
			lwsz = raw_lwsz;
			break;
		end

	end

	ind = find(res_raw(:, 1) < lwsz);
	res_raw(ind, :) = [];

end

if nsp_flag

	while 1

		raw_sz_cf = input('please type in lower cutoff of liposome size for removing non-specific detection: ');

		if raw_sz_cf <= 0 || ~isnumeric(raw_sz_cf)
			disp('input value should be a positive number, please retype.');
			continue;
		else
			sz_cf = raw_sz_cf;
			break;
		end

	end

	while 1

		raw_int_cf = input('please type in density cutoff for removing non-specific detection: ');

		if raw_int_cf <= 0 || ~isnumeric(raw_int_cf)
			disp('input value should be a positive number, please retype.');
			continue;
		else
			int_cf = raw_int_cf;
			break;
		end

	end

	ind2 = find(res_raw(:, 1) < sz_cf & res_raw(:, 2) < int_cf);
	res_raw(ind2, :) = [];

end

out = res_raw;





