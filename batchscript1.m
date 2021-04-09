clear;

list_pep = dir(strcat('E:\data\fig.s for paper\SLiC\181209 SLiC\PD1SDP2T_K\1\pep\*.tif'));
list_pep_num = length(list_pep);

list_lip = dir(strcat('E:\data\fig.s for paper\SLiC\181209 SLiC\PD1SDP2T_K\1\lip\*.tif'));
list_lip_num = length(list_lip);

aa = {};
counter = [];

for j = 1:list_pep_num
	img_np = list_pep(j).name;
    cd('E:\data\fig.s for paper\SLiC\181209 SLiC\PD1SDP2T_K\1\pep\');
	g_pep = double(imread(strcat(img_np)));
	img_nl = list_lip(j).name;
    cd('E:\data\fig.s for paper\SLiC\181209 SLiC\PD1SDP2T_K\1\lip\');
	g_lip = double(imread(strcat(img_nl)));

	g_pep_bp = bpfilter(g_pep, 1, 5, 0.5);
	det_pep = spotmulsz_ver2forbatch(g_pep_bp, 1);

	if ~isempty(det_pep)

		out = spotscorr_ver2forbatch(g_pep_bp, g_lip, det_pep, 4.97, 5.65);

	else
		out = [];
	end

	if ~isempty(out)

		out = sliccutoff_nonint(out, 1, 1);

	end

	aa{j} = out;

	counter(j) = length(aa{j});

	cd('E:\data\fig.s for paper\SLiC\181209 SLiC\PD1SDP2T_K\1\');

	if ~isempty(out);

		xlswrite(['res', num2str(j), '.xls'], aa{j});

	end
end

xlswrite('counts.xls', counter');

disp('Done.');


