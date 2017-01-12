function [] = run_algorithm_comparison()
% [] = run_algorithm_comparison()
%
% Script for comparing TSP optimizers:
% generates convergence plots and TSP solution plots
%
% Naming convention (optimizers then automatically detected):
% 	lastname1_lastname2_sa.m lastname_ga.m
%	  lastname_sa() lastname_ga() 	
%
% Author: Hao Wang, Koen van der Blom

	% Comparison setup
	runs_per_optimizer = 20; 
	eval_budget = 15000;
    
    % problem setup: [order, objfunc]
	problems = { ...
        12, 'eval_square'; ...
        9, 'eval_semi_cube'; ...
        9, 'eval_cube'
    };
    n_problems = size(problems, 1);

	% Initialize statistics  
	statistics([], eval_budget);  

	% Load optimizers
	sa_optimizers = dir('*sa.m');
	sa_optimizers = {sa_optimizers.name};
	ga_optimizers = dir('*ga.m');
	ga_optimizers = {ga_optimizers.name};
    
	optimizers = [sa_optimizers ga_optimizers];
	optimizers = strrep(optimizers, '.m', '');
	optimizer_names = strrep(optimizers, '_', '\_');
    
	if isempty(optimizers)
		disp('No optimizers detected!');
		return
  	else
		fprintf('%s', [num2str(length(optimizers)), ' optimizers detected: '])
		disp(optimizers);
	end

	% Output settings
	colors = {'b', 'r', 'g', 'm', 'c', 'k', '-.b', '-.r', '-.g', '-.m', '-.c', '-.k'};
	tab = '    ';
    myStyle = hgexport('factorystyle');
    myStyle.Format = 'jpeg';
    myStyle.Resolution = 300;
	results_dir = 'Results/';
    
	if (~exist(results_dir, 'dir')); mkdir(results_dir); end

	% For each problem run all optimizers, generate convergence plot, and plot best square (cube) found per optimizer
    hist_best_so_far = zeros(n_problems, length(optimizers), length(runs_per_optimizer), eval_budget);
    elapsed = zeros(n_problems, length(optimizers), length(runs_per_optimizer));
	for i = 1 : n_problems
        
        [order, prob_str] = problems{i, :};
        prob_name = strrep(prob_str, 'eval_', '');
        if isempty(regexp(prob_name, 'cube', 'once'))
            dim = 2;
        else
            dim = 3;
        end
        
		disp(['Test problem ' num2str(i), '/', num2str(n_problems), ' (',  prob_name, ')'])
		hist_x_opt = zeros(length(optimizers), length(runs_per_optimizer), order^dim);
        
		for j = 1 : length(optimizers)
			disp([tab, 'Optimizer ', num2str(j), '/', num2str(length(optimizers)), ' (', optimizers{j}, ') on (',  prob_name, '):'])
			optimizer = str2func(optimizers{j});
            
			for k = 1 : runs_per_optimizer
				run_file = [results_dir, prob_name, '_', cell2mat(optimizers(j)), '_', num2str(k), '.mat'];
                
% obsolete....                
% 				if (exist(run_file,'file'))
% 					fprintf('%s', [tab, tab, 'Loading file ''', run_file, ''': ']);
% 					load(run_file);
% 				else
% 					fprintf('%s', [tab, tab, 'Executing run ', num2str(k), '/', num2str(runs_per_optimizer), ': ']);
% 					tic;
% 					[stat.xopt, stat.fopt] = optimizer(dim, eval_budget, str2func(prob_str));
% 					stat.elapsed = toc;
% 					stat.hist_best_so_far = statistics([], eval_budget);
% 					save(run_file, 'stat');
%                 end

                fprintf('%s', [tab, tab, 'Executing run ', num2str(k), '/', num2str(runs_per_optimizer), ': ']);
                tic;
                [stat.xopt, stat.fopt, stat.hist_best_so_far] = optimizer(dim, eval_budget, str2func(prob_str));
%                 [stat.xopt, stat.fopt] = optimizer(dim, eval_budget, str2func(prob_str));

                stat.elapsed = toc;
%                 stat.hist_best_so_far = statistics([], eval_budget);
                
                save(run_file, 'stat');
                
				hist_x_opt(j, k, :) = stat.xopt;
				hist_best_so_far(i, j, k, :) = stat.hist_best_so_far(1:eval_budget);
				elapsed(i, j, k) = stat.elapsed;
				fprintf('fopt=%f, elapsed=%f\n', hist_best_so_far(i, j, k, eval_budget), stat.elapsed)
                
            end
            
			if runs_per_optimizer > 1
 				fprintf('%smedian: fopt=%f, elapsed=%f\n', [tab, tab], median(hist_best_so_far(i, j,: , eval_budget)), median(elapsed(i, j, :)));
			end
		end

		% ---------------------------------------------------------------------------
		% Convergence plot
		% ---------------------------------------------------------------------------
		fig = figure;
		for j = 1:length(optimizers)
			plot_hist_best_so_far = squeeze(hist_best_so_far(i,j,:,1:eval_budget));
			if runs_per_optimizer > 1
				plot_hist_best_so_far = median(plot_hist_best_so_far);
			end
			plot([1:eval_budget], plot_hist_best_so_far, char(colors(j)),'LineWidth', 1.5)
			hold on
		end
		legend(optimizer_names, 'Location', 'Best');
		grid on
		ylabel('fitness', 'FontWeight', 'Bold', 'FontSize', 10);
		xlabel('evaluations', 'FontWeight', 'Bold', 'FontSize', 10);
		set(gca, 'FontWeight', 'Bold', 'FontSize', 10);
		title([prob_name, ': convergence'], 'FontWeight', 'Bold', 'fontsize', 12);
		savefile_plot = [results_dir, prob_name, '_convergence'];
		hgexport(fig, savefile_plot, myStyle);

		% ---------------------------------------------------------------------------
		% Plot the best magic square (cube) found per optimizer
		% --------------------------------------------------------------------------
		for j = 1:length(optimizers)
			fig = figure('Position', [100, 100, 1000, 1000]);
			[best_score, best_index] = min(hist_best_so_far(i,j,:,end));
            if dim == 2
                plot_square(hist_x_opt(j, best_index, :));
            elseif dim == 3
                plot_cube(hist_x_opt(j, best_index, :));
            end
            
            % obsolete...
			% set(gca, 'FontWeight', 'Bold', 'FontSize', 10);
			% title([cell2mat(optimizer_names(j)), ' : ', prob_name, ' (', num2str(best_score), ')'], 'FontWeight', 'Bold', 'fontsize', 12);
			savefile_plot = [results_dir, prob_name, '_', cell2mat(optimizers(j)), '_best_ever_found'];
            hgexport(fig, savefile_plot, myStyle);
		end

		close all

	end

end


function [] = plot_square(x)

    n = sqrt(length(x));
    X = reshape(x, n, n);
    t = uitable('Data', X, 'ColumnWidth', repmat({40}, 1, n), 'Position', [191.5, 324.5, 517, 251]);
    tmp = get(t, 'Position');
    tmp = [tmp(1:2) 0 0] + get(t, 'Extent');   
    set(t, 'Position', tmp);

end



function [] = plot_cube(x)

    n = nthroot(length(x), 3);
    
    X = reshape(x, [n, n, n]);
    
    h_space = (1000 - 304 * 3) / 4;
    v_space = (1000 - 194 * 3) / 4;
    width = 304;
    height = 194;
    
    for i = 1 : n
        c = mod(i, 3);
        if c == 0; c = 3; end
        r = ceil(i / 3);
        left = h_space * c + width * (c-1);
        bottom = v_space * (3 - r + 1) + height * (3 - r);
        t = uitable('Data', X(:, :, i), 'ColumnWidth', repmat({30}, 1, n), 'Position', [left, bottom, width, height]);
        tmp = get(t, 'Position');
        tmp = [tmp(1:2) 0 0] + get(t, 'Extent');   
        set(t, 'Position', tmp);
    end

end
