function [C, T, converged, clust, Cf, varargout] = driver_acdc(Tmax, varargin)
%function [C, T, converged, clust, Cf, labels_ch, clust_flux, J_out, flux, outflux_matrix, out, sources] = driver_acdc(Tmax, {C0, {T0}}, {'File', Filename, {...}}))
%function [C, T, converged, clust, Cf, ~, ~, J_out] = driver_acdc(Tmax, {C0, {T0}}, 'no_dofluxes',  {'Option', {Value}, {...}}))
%function [C, T, converged, clust, Cf] = driver_acdc(Tmax, {C0, {T0}}, 'no_fluxes',  {'Option', {Value}, {...}}))
% 
% Input parameters:
% Tmax:		     duration of the simulation in seconds (added to T0), or vector containing output time points (not added to T0).
% C0:			 an optional row vector (or matrix) containing the initial concentrations (on the last row).
% T0:			 if C0 contains concentrations as a function of time, the corresponding times should be given in the column vector T0.
% 'No_dofluxes':  the driver does not call dofluxes.m and no related output is printed
% 'No_fluxes':    the driver does not call dofluxes.m or formationrate.m and no related output is printed
% 'Sources_in', filename:    input file containing source terms, concentrations to be kept constant and initial concentrations
% 'Sources_out', filename:   output file for monomer source terms, only used if steady-state is reached
% 'Constants_out', filename: output file for constant monomer concentrations, only used if steady-state is reached
% 'C_neutr', filename:	    output file for concentrations of neutrals
% 'C_neg', filename:		    output file for concentrations of negatives
% 'C_pos', filename:		    output file for concentrations of positives
% 'Fluxes', filename:	    output file for fluxes
% 'Outmat', filename:	    output file for matrix containing collisions going out from the system
% 'Cluster_data', filename:	    output mat-file for cluster names, masses, diameters etc.
% 'Cfun', {cluster function variable1 variable2 ...}: function giving the concentration of some molecule or cluster as a function of time or other concentrations, example: 'Cfun', {'1A' 'acidfunction' 'T'} (acidfunction.m gives acid concentration as a function of time) or 'Cfun', {'1D' 'dma_from_ammonia' '1N' 'T'} (dma_from_ammonia.m gives DMA concentration as a function of ammonia concentration and time)
% 
% Output parameters:
% C:			  cluster distribution in cm^(-3) as a function of time T
% T:			  time points (seconds) corresponding to C
% converged:	  1 if the simulation has converged, and 0 if not, and -1 if there are negative concentrations lower than the threshold value -1e-12 cm^(-3)
% clust:		  cluster names
% Cf:          final cluster distribution in cm^(-3)
% labels_ch:	  names of the neutral and charged forms of all molecules, removed or extra protons and generic ions
% clust_flux:	  names for all fluxes, including also non-clusters (sources, losses etc.)
% J_out:       particle formation rate (rate out of the system) in cm^(-3) s^(-1) as a function of time
% flux:        net fluxes from i to j (cm^(-3) s^(-1)) are given in the matrix elements flux(i,j)
%              in collisions/fissions involving two similar clusters, the flux is from the point of view of the
%              SMALLER cluster, i.e. must be divided by 2 when considering the point of view of the larger cluster
% outflux_matrix:	  fluxes going out of the system by collision of clusters i and j (cm^(-3) s^(-1)) 
%     		  are given in the matrix elements outflux_matrix(i,j)
%     		  for i=j, the flux is from the point of view of the outgrowing clusters (i.e. does NOT need to be divided by 2)
% out:         fluxes from i out to a cluster with j mol1 (cm^(-3) s^(-1)) are in out(i,j+1,k) with
%     		  k=1: neutral, k=2: neg, k=3: pos, k=4: recombination
%     		  in collisions involving two similar clusters, the flux is from the point of view of 
%     		  the outgrowing clusters (i.e. does NOT need to be divided by 2)
% sources:	  sources of each species (cm^(-3) s^(-1))


% Cluster definitions etc.

	% molecule names
	labels = {'A' 'B' 'N' 'P'};
	% neutral, negative and positive species
	labels_ch = {'A' 'B' '' ; 'N' '' '1N1P' ; '' 'B-1A' 'P' ; '' 'neg' 'pos'};
	% cluster names
	clust = {'1A' '2A' '3A' '4A' '5A' '1N' '1A1N' '2A1N' '3A1N' '4A1N' '5A1N' '2N' '1A2N' '2A2N' '3A2N' '4A2N' '5A2N' '3N' '1A3N' '2A3N' '3A3N' '4A3N' '5A3N' '4N' '1A4N' '2A4N' '3A4N' '4A4N' '5A4N' '4A5N' '5A5N' '1B' '1A1B' '2A1B' '3A1B' '4A1B' '1B1N' '1A1B1N' '2A1B1N' '3A1B1N' '4A1B1N' '1B2N' '1A1B2N' '2A1B2N' '3A1B2N' '4A1B2N' '3A1B3N' '4A1B3N' '4A1B4N' '1N1P' '1A1N1P' '2A1N1P' '2N1P' '1A2N1P' '2A2N1P' '3A2N1P' '3N1P' '1A3N1P' '2A3N1P' '3A3N1P' '4A3N1P' '2A4N1P' '3A4N1P' '4A4N1P' '3A5N1P' '4A5N1P' '5A5N1P' 'neg' 'pos'};
	% name vector containing also the names of fluxes other than clusters (sources, losses etc.)
	clust_flux = [clust 'source' 'coag' 'wall' 'dil' 'rec' 'out_neu' 'out_neg' 'out_pos' 'bound'];
	charging_states = [1 1 1 1 1 1 1 1 1 1 ...
		1 1 1 1 1 1 1 1 1 1 ...
		1 1 1 1 1 1 1 1 1 1 ...
		1 2 2 2 2 2 2 2 2 2 ...
		2 2 2 2 2 2 2 2 2 3 ...
		3 3 3 3 3 3 3 3 3 3 ...
		3 3 3 3 3 3 3 2 3]; % charging states (1 = neutral, 2 = negative, 3 = positive)
	% dry values
	diameters = [0.55 0.70 0.80 0.88 0.95 0.43 0.63 0.75 0.84 0.91 ...
		0.98 0.54 0.69 0.79 0.87 0.94 1.00 0.62 0.74 0.83 ...
		0.91 0.97 1.03 0.68 0.78 0.87 0.94 1.00 1.05 1.02 ...
		1.07 0.55 0.70 0.80 0.88 0.95 0.63 0.75 0.84 0.91 ...
		0.97 0.69 0.79 0.87 0.94 1.00 0.97 1.03 1.05 0.43 ...
		0.63 0.75 0.54 0.69 0.79 0.87 0.62 0.74 0.83 0.91 ...
		0.97 0.87 0.94 1.00 0.96 1.02 1.07 0.45 0.39];
	% dry values
	mobility_diameters = [0.85 1.00 1.10 1.18 1.25 0.73 0.93 1.05 1.14 1.21 ...
		1.28 0.84 0.99 1.09 1.17 1.24 1.30 0.92 1.04 1.13 ...
		1.21 1.27 1.33 0.98 1.08 1.17 1.24 1.30 1.35 1.32 ...
		1.37 0.85 1.00 1.10 1.18 1.25 0.93 1.05 1.14 1.21 ...
		1.27 0.99 1.09 1.17 1.24 1.30 1.27 1.33 1.35 0.73 ...
		0.93 1.05 0.84 0.99 1.09 1.17 0.92 1.04 1.13 1.21 ...
		1.27 1.17 1.24 1.30 1.26 1.32 1.37 0.75 0.69]; % mobility diameter (2*r+0.3nm)*sqrt(1+28.8u/m)
	masses = [98.08 196.16 294.24 392.32 490.40 17.04 115.12 213.20 311.28 409.36 ...
		507.44 34.08 132.16 230.24 328.32 426.40 524.48 51.12 149.20 247.28 ...
		345.36 443.44 541.52 68.16 166.24 264.32 362.40 460.48 558.56 477.52 ...
		575.60 97.08 195.16 293.24 391.32 489.40 114.12 212.20 310.28 408.36 ...
		506.44 131.16 229.24 327.32 425.40 523.48 442.44 540.52 557.56 18.04 ...
		116.12 214.20 35.08 133.16 231.24 329.32 52.12 150.20 248.28 346.36 ...
		444.44 265.32 363.40 461.48 380.44 478.52 576.60 32.00 19.02]; % mass in amus
	monomers = [1 6 32 50];
	nonmonomers = [2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69];
	neutrals = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
		21 22 23 24 25 26 27 28 29 30 31];
	negatives = [32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 ...
		47 48 49 68];
	positives = [50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 ...
		65 66 67 69];
	neutral_monomers = [1 6];
	negative_monomers = [32];
	positive_monomers = [50];
	neutral_clusters = [2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
		21 22 23 24 25 26 27 28 29 30 31];
	negative_clusters = [33 34 35 36 37 38 39 40 41 42 43 44 45 46 ...
		47 48 49];
	positive_clusters = [51 52 53 54 55 56 57 58 59 60 61 62 63 64 ...
		65 66 67];
	A_in_clust = [1 2 3 4 5 0 1 2 3 4 ...
		5 0 1 2 3 4 5 0 1 2 ...
		3 4 5 0 1 2 3 4 5 4 ...
		5 1 2 3 4 5 1 2 3 4 ...
		5 1 2 3 4 5 4 5 5 0 ...
		1 2 0 1 2 3 0 1 2 3 ...
		4 2 3 4 3 4 5 0 0];
	table_0 = {[2 1] [3 1] [4 1] [5 1] [6 1] [1 2] [2 2] [3 2] [4 2] [5 2] [6 2] [1 3] [2 3] [3 3] [4 3] [5 3] [6 3] [1 4] [2 4] [3 4] [4 4] [5 4] [6 4] [1 5] [2 5] [3 5] [4 5] [5 5] [6 5] [5 6] [6 6]};
	table_neg = {[2 1] [3 1] [4 1] [5 1] [6 1] [2 2] [3 2] [4 2] [5 2] [6 2] [2 3] [3 3] [4 3] [5 3] [6 3] [5 4] [6 4] [6 5]};
	table_pos = {[1 2] [2 2] [3 2] [1 3] [2 3] [3 3] [4 3] [1 4] [2 4] [3 4] [4 4] [5 4] [3 5] [4 5] [5 5] [4 6] [5 6] [6 6]};
	% criteria for outgrowing neutrals
	nmols_out_neutral = reshape([6 0 5 0],[1 4]);
	% criteria for outgrowing negatives
	nmols_out_negative = reshape([5 1 3 0],[1 4]);
	% criteria for outgrowing positives
	nmols_out_positive = reshape([5 0 6 1],[1 4]);


% Initializing concentrations etc.

	C0 = zeros(1,78);

	isconst = zeros(69,1);
	isfitted = zeros(69,1);
	C_fit = cell(69,2);
	source = zeros(69,1);

	J_out = nan;
	flux = nan;
	outflux_matrix = nan;
	out = nan;
	sources = nan;

% Checking what input we have and processing it

	if size(Tmax,1)>size(Tmax,2), Tmax = Tmax'; end
	fixed_time = 1;
	if nargout<=7
		no_fluxes = 1;
		no_j = 1;
		calc_sources = 0;
	elseif nargout<=8
		no_fluxes = 1;
		no_j = 0;
		calc_sources = 0;
	else
		no_fluxes = 0;
		no_j = 0;
		calc_sources = 1;
	end
	sources_in = 0;
	sources_out = 0;
	consts_out = 0;
	opt_str = '';
	c_neutr = 0;
	c_neg = 0;
	c_pos = 0;
	fl_out = 0;
	outmat = 0;
	cluster_data = 0;
	C_fun_str = {};
	n_fun = [];
	C0_in = 0;
	CSfactor = 1;
	WLfactor = 1;
	if ~isempty(varargin)
		i = 0;
		if isnumeric(varargin{1})	% Initial concentrations cm^3 -> m^3
			C00 = varargin{1}*1e6;
			C0_in = 1;
			i = 1;
			if size(C00,2)==69
				C00 = [C00 zeros(size(C00,1),9)];
				C0 = C00(end,:);
				if (length(varargin)>1 && isnumeric(varargin{2}))	% Time (vector) corresponding to C00
					T00 = varargin{2};
					if (isvector(T00) && numel(T00)==size(C00,1))
						i = 2;
						if size(T00,2)>size(T00,1), T00 = T00'; end
						if numel(Tmax)==1
							Tmax = Tmax+T00(end);
						end;
					else
						error('Sizes of initial concentration array and corresponding time array don''t match.')
					end
				else
					T00 = 0;
				end
			else
				error('Bad initial concentrations.')
			end
		else
			T00 = 0;
			C00 = C0;
		end
		% Check for single keywords in the input arguments
		if ismember(1,strcmpi('fixed_time',varargin))
			fixed_time = 1;
			varargin = varargin(~strcmpi('fixed_time',varargin));
		elseif ismember(1,strcmpi('repeat',varargin))
			fixed_time = 0;
			varargin = varargin(~strcmpi('repeat',varargin));
		elseif ismember(1,strcmpi('no_dofluxes',varargin))
			no_fluxes = 1;
			varargin = varargin(~strcmpi('no_dofluxes',varargin));
		elseif ismember(1,strcmpi('no_fluxes',varargin))
			no_fluxes = 1;
			no_j = 1;
			varargin = varargin(~strcmpi('no_fluxes',varargin));
		end
		% Go through the paired input arguments
		for j = i+1:2:length(varargin)
			if strcmpi(varargin{j},'Sources_in')
				sources_in = 1;
				fname1 = varargin{j+1};
			elseif strcmpi(varargin{j},'Sources_out')
				sources_out = 1;
				fname2 = varargin{j+1};
				no_fluxes = 0;
				calc_sources = 1;
			elseif strcmpi(varargin{j},'Constants_out')
				consts_out = 1;
				fname2c = varargin{j+1};
			elseif strcmpi(varargin{j},'Options')
				opt_str = varargin{j+1};
			elseif strcmpi(varargin{j},'C_neutr')
				c_neutr = 1;
				fname3 = varargin{j+1};
			elseif strcmpi(varargin{j},'C_neg')
				c_neg = 1;
				fname4 = varargin{j+1};
			elseif strcmpi(varargin{j},'C_pos')
				c_pos = 1;
				fname5 = varargin{j+1};
			elseif strcmpi(varargin{j},'Fluxes')
				fl_out = 1;
				fname6 = varargin{j+1};
				no_fluxes = 0;
				calc_sources = 1;
			elseif strcmpi(varargin{j},'CSfactor')
				CSfactor = varargin{j+1};
			elseif strcmpi(varargin{j},'WLfactor')
				WLfactor = varargin{j+1};
			elseif strcmpi(varargin{j},'Outmat')
				outmat = 1;
				fname7 = varargin{j+1};
				no_fluxes = 0;
			elseif strcmpi(varargin{j},'Cluster_data')
				cluster_data = 1;
				fname8 = varargin{j+1};
				no_fluxes = 0;
				calc_sources = 1;
			elseif strcmpi(varargin{j},'Cfun')
				n_fun = [n_fun find(strcmp(clust,varargin{j+1}{1}))];
				isconst(n_fun(end)) = 1;
				if length(varargin{j+1})<2 || isempty(varargin{j+1}{2})
					error(['The function related to cluster ',varargin{j+1}{1},' was not given.'])
				end
				C_fun_str{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];
				C_fun_str_2{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];
				C_fun_str_3{n_fun(end)} = ['1e6*', varargin{j+1}{2},'('];
				for k = 3:length(varargin{j+1})
					if strcmpi(varargin{j+1}{k},'t')
						C_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, 't'];
						C_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, 'T(j)'];
						C_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, 'T'];
					else
						C_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, '1e-6*c(', num2str(find(strcmp(clust,varargin{j+1}{k}))), ')'];
						C_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, '1e-6*C(j,', num2str(find(strcmp(clust,varargin{j+1}{k}))), ')'];
						C_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, '1e-6*C(j,', num2str(find(strcmp(clust,varargin{j+1}{k}))), ')'];
					end
					if k<length(varargin{j+1})
						C_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, ','];
						C_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, ','];
						C_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, ','];
					end
				end
				C_fun_str{n_fun(end)} = [C_fun_str{n_fun(end)}, ')'];
				C_fun_str_2{n_fun(end)} = [C_fun_str_2{n_fun(end)}, ')'];
				C_fun_str_3{n_fun(end)} = [C_fun_str_3{n_fun(end)}, ')'];
			end
		end
	end
	if ~isempty(n_fun)
		n_fun = sort(n_fun);
		dC_fun_str = 'dcdt = @(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs) equations_acdc(t,[c(1:';
		for i = 1:length(n_fun)
			dC_fun_str = [dC_fun_str, num2str(n_fun(i)-1),');',C_fun_str{n_fun(i)}, ';c(', num2str(n_fun(i)+1), ':'];
		end
		dC_fun_str = [dC_fun_str, 'end)],K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs);'];
		eval(dC_fun_str);
	else
		dcdt = @(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs) equations_acdc(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs);
	end

	if numel(Tmax)>1
		if Tmax(1)==T00(end)
			Tmax = Tmax(2:end);
		end
	end

	if cluster_data
		save(fname8, 'clust', 'diameters', 'mobility_diameters', 'masses', 'monomers', 'nonmonomers', 'neutrals', 'negatives', 'positives', 'clust_flux','A_in_clust')
	end

% If input source file was given, reading its content and converting cm^3 -> m^3
	if sources_in

		sfid = fopen(fname1, 'r');
		if sfid == -1
		error('Source file not found!')
		end
		A = textscan(sfid, '%s %s %f %s', 'commentStyle', '%');
		fclose(sfid);
		A2 = A{2}(cellfun('isempty',regexp(A{1},'^(wall|coag)')));
		if any(~ismember(A2,clust)), warning(['Discarding erroneous cluster names in ',fname1]), end

		for i=1:69
			j = find(strcmp(A{2},clust{i}));
			defined = 0;
			for k=1:length(j)
				if (~C0_in && strncmp(A{1}(j(k)),'i',1))
					if (defined == 1 || defined == 3), error('Source term multiply defined.'), end
					C0(i) = A{3}(j(k))*1e6;
					defined = 1;
				elseif strncmp(A{1}(j(k)),'s',1)
					if (defined > 1), error('Source term multiply defined.'), end
					source(i) = A{3}(j(k))*1e6;
					defined = 2;
				elseif strncmp(A{1}(j(k)),'c',1)
					if (defined > 0), error('Source term multiply defined.'), end
					if length(A{4})<j(k) || isempty(A{4}{j(k)})
						isconst(i) = 1;
						C0(i) = A{3}(j(k))*1e6;
					else
						C_fit{i,1} = A{3}(j(k))*1e6;
						C_fit{i,2} = [];
						A4 = regexp(A{4}{j(k)},'-','split');
						A4(cellfun('isempty',A4)) = [];
						if any(~ismember(A4,clust)), warning(['Discarding erroneous cluster names in ',fname1]), end
						for i2=1:69
							for l=1:length(find(strcmp(A4,clust{i2})))
								C_fit{i,2} = [C_fit{i,2},i2];
							end
						end
						if ~isempty(C_fit{i,2})
							isfitted(i) = 1;
						else
							isconst(i) = 1;
						end
						C0(i) = A{3}(j(k))*1e6-sum(C0(C_fit{i,2}));
					end
					defined = 3;
				end
			end
		end

		C00(end,:) = C0;
		j = find(strncmp(A{1},'wall',4));
		if length(j) > 1
			error('Wall loss enhancement multiply defined.')
		end
		if ~isempty(j),
			fwl_in = A{3}(j);
		end
		j = find(strncmp(A{1},'coag',4));
		if length(j) > 1
			error('Coagulation loss enhancement multiply defined.')
		end
		if ~isempty(j),
			fcs_in = A{3}(j);
		end
	elseif C0_in == 0
		error('Input concentrations or source terms must be given either in the vector C0 or the file ''Sources_in''!');
	end

% Setting the options for the solver
	eval(['options = odeset(',opt_str,');']);

% Reading in the collision and evaporation rates, coagulation losses and wall losses (enhancement factors for ion losses included in the loss rate files)
	get_coll;
	get_evap;
	get_cs;
	get_wl;

	CS = CS*CSfactor;
	WL = WL*WLfactor;
	if exist('fcs_in','var')
		fcs = fcs_in;
	end
	if exist('fwl_in','var')
		fwl = fwl_in;
	end

% Solving the birth-death equations
	converged = 0;
	i = 0;
	Tfin = Tmax;
	if fixed_time
		imax = 1;
	else
		imax = 5;
		Tadd = Tmax-T00(end);
	end
	while ~converged
		i = i+1;
% Checking for negative concentrations
		if min(min(C00))<-1e-6
			converged = -1;
			Cf = nan(size(C0));
			C = C00*1e-6;
			T = T00;
			J_out = nan;
			flux = nan;
			outflux_matrix = nan;
			out = nan;
			sources = nan;
			return;
		end
		if i>imax
% Not converged but exiting anyway
			if nargout < 3
				disp('Not converging, try a larger Tmax');
			end
			break;
		elseif i>1
% Adding more time in order to reach convergence
			Tfin = T00(end)+Tadd;
		end
		[T,C] = ode15s(@(t,c) dcdt(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs), [T00(end) Tfin], C0, options);
		if i == 1 && min(min(C)) < 0;
			[T5,C5] = ode15s(@(t,c) dcdt(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs), [T00(end) T00(end)+(Tfin(1)-T00(end))*1e-6], C0, options);
			[T,C] = ode15s(@(t,c) dcdt(t,c,K,E,WL,CS,source,isconst,isfitted,C_fit,fwl,fcs), [T5(end) Tfin], C5(end,:), options);
			if numel(Tmax)==1
				T = [T5; T(2:end)];
				C = [C5; C(2:end,:)];
			end
		end
		Cf = C(end,:);
		if max(abs(Cf(1:69)-C0(1:69))./max(C0(1:69),1e-6))<1e-6
			converged = 1;
		elseif length(T)>2 && max(abs(Cf(1:69)-C(floor(length(T)*.75),1:69))./max(C(floor(length(T)*.75),1:69),1e-6))<1e-6
			converged = 1;
		end
		if any(isfitted)
			for j = 1:69
				if isfitted(j)
					if abs(C_fit{j,1}-Cf(j)-sum(Cf(C_fit{j,2})))/C_fit{j,1} > 1e-6
						Cf([j C_fit{j,2}]) = Cf([j C_fit{j,2}])*C_fit{j,1}/sum(Cf([j C_fit{j,2}]));
						converged = 0;
					end
				end
			end
		end
		C0 = Cf;
		C00 = [C00; C(2:end,:)];
		T00 = [T00; T(2:end)];
	end
	if ~C0_in && numel(Tmax)==1
		C = C00(2:end,:);
		T = T00(2:end);
	else
		C = C00;
		T = T00;
	end
	for i=n_fun
		C(:,i) = eval(C_fun_str_3{i});
	end

% Printing table of neutral concentrations
	if(c_neutr)
		t0 = zeros(5,6);
		for i = 1:length(neutrals)
			t0(table_0{i}(1),table_0{i}(2)) = neutrals(i);
		end
		cfid = fopen(fname3, 'w');
		fprintf(cfid,'%% %-9s %-9s %-9s %-9s %-9s %-9s \n', '0N0D', '1N0D', '2N0D', '3N0D', '4N0D', '5N0D');
		for i=1:size(t0,1)
			for j=1:size(t0,2)
				if t0(i,j)<1
					fprintf(cfid,'	  nan ');
				else
					fprintf(cfid,'%.3e ',Cf(t0(i,j))*1e-6);
				end
			end
			fprintf(cfid,'\n');
		end
		fclose(cfid);
	end

% Printing table of negative concentrations
	if(c_neg)
		tn = zeros(4,5);
		for i = 1:length(negatives)
			tn(table_neg{i}(1),table_neg{i}(2)) = negatives(i);
		end
		cfid = fopen(fname4, 'w');
		fprintf(cfid,'%% %-9s %-9s %-9s %-9s %-9s \n', '0N0D', '1N0D', '2N0D', '3N0D', '4N0D');
		for i=1:size(tn,1)
			for j=1:size(tn,2)
				if tn(i,j)<1
					fprintf(cfid,'	  nan ');
				else
					fprintf(cfid,'%.3e ',Cf(tn(i,j))*1e-6);
				end
			end
			fprintf(cfid,'\n');
		end
		fclose(cfid);
	end

% Printing table of positive concentrations
	if(c_pos)
		tp = zeros(5,6);
		for i = 1:length(positives)
			tp(table_pos{i}(1),table_pos{i}(2)) = positives(i);
		end
		cfid = fopen(fname5, 'w');
		fprintf(cfid,'%% %-9s %-9s %-9s %-9s %-9s %-9s \n', '0N0D', '1N0D', '2N0D', '3N0D', '4N0D', '5N0D');
		for i=1:size(tp,1)
			for j=1:size(tp,2)
				if tp(i,j)<1
					fprintf(cfid,'	  nan ');
				else
					fprintf(cfid,'%.3e ',Cf(tp(i,j))*1e-6);
				end
			end
			fprintf(cfid,'\n');
		end
		fclose(cfid);
	end

	if ~no_j
		J_out = nan(1,length(T));
		for i=1:length(T)
			J_out(i) = formationrate(C(i,:)*1e-6,K);
		end
	end
	if ~no_fluxes
		[flux, ~, sources, ~, flux_out, ~] = dofluxes(C(end,:),K,E,WL,CS,source,calc_sources,fwl,fcs);
		% Back to cubic centimeters
		flux = flux*1e-6;
		sources = sources*1e-6;
		flux_out = flux_out*1e-6;
	end

	C = C*1e-6;
	Cf = Cf*1e-6;
	% printing out source file for constant monomer concentrations
	if consts_out && converged
		cfid = fopen(fname2c, 'w');
		for i=monomers
			fprintf(cfid,'constant %s %f\n',clust{i}, Cf(i));
		end
		fclose(cfid);
	end

	if ~no_fluxes
		% printing out source file for monomer sources
		if sources_out && converged
			sfid = fopen(fname2, 'w');
			for i=monomers
				fprintf(sfid,'source %s %f\n',clust{i}, sources(i));
			end
			fclose(sfid);
		end

		out = zeros(69,13,4);
		outflux_matrix = zeros(69,69);
		for i=1:69
			for j=1:69
				k = A_in_clust(i)+A_in_clust(j)+1;
				out(i,k,:) = out(i,k,:)+flux_out(i,j,:);
				outflux_matrix(i,j) = sum(flux_out(i,j,:));
			end
		end
	% Printing the flux matrix
		if(fl_out)
			fid=fopen(fname6,'w');
			fprintf(fid,'%% ');
			for i=1:length(clust_flux)
				fprintf(fid,'%s ',clust_flux{i});
			end
			fprintf(fid,'\n');
			for i=1:size(flux,1)
				for j=1:size(flux,2)
					fprintf(fid,'%.6e ',flux(i,j));
				end
				fprintf(fid,'\n');
			end
			fclose(fid);
		end

	% Printing the outflux matrix
		if(outmat)
			fid=fopen(fname7,'w');
			fprintf(fid,'%% ');
			for i=1:length(clust)
				fprintf(fid,'%s ',clust{i});
			end
			fprintf(fid,'\n');
			for i=1:size(outflux_matrix,1)
				for j=1:size(outflux_matrix,2)
					fprintf(fid,'%.6e ',outflux_matrix(i,j));
				end
				fprintf(fid,'\n');
			end
			fclose(fid);
		end

	end
	C = C(:,1:69);
	Cf = Cf(1:69);
	if nargout>5
		output_opt={labels_ch, clust_flux, J_out, flux, outflux_matrix, out, sources};
		varargout=output_opt(1:nargout-5);
	end
end