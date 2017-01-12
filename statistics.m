function [hist_best_so_far] = statistics(f, eval_budget)
% Retrieve and reset:	statistics([], eval_budget)
% Update stats: 	statistics(f)

	persistent p_eval_count p_eval_budget p_hist_best_so_far

	% Retrieve and reset
	if isempty(f)
		p_eval_count = 0;
		p_eval_budget = eval_budget;

		hist_best_so_far = p_hist_best_so_far;
	% Update stats
	else
		if p_eval_count == 0
			p_hist_best_so_far = zeros(p_eval_budget, 1);
			f_best = f;
		else 
			f_best = p_hist_best_so_far(p_eval_count);
			if f < f_best % M I N I M I Z A T I O N 
				f_best = f;
			end
		end

		p_eval_count = min(p_eval_budget, p_eval_count + 1);
		p_hist_best_so_far(p_eval_count:p_eval_budget) = f_best;
	end
end
