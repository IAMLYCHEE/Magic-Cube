function [xopt_pheno, fopt,hist_fitness] = yicheng_li_ga(dim, eval_budget, fitness_func)

%   skeleton for simulated anealing algorithm. For implementation, we 
%   recommend to use a vector of integers between [1, n^2] (for magic 
%   square) as the representation (adapt such a representation to the 
%   magic cube problem by yourself). In addition, you have to come up 
%   with the mutation operator by yourself :)
%
%   Parameters
%   ----------
%       dim : problem dimension, which should be either 2 (square) or 3 (cube)
%       eval_budget : int, the function evaluation budget
%       fitness_func : function handle, you should use one of the
%           evaluation function provided
%           
%   Output
%   ----------
%       xopt : array, the final solution vector found by the algorithm
%       fopt : double, the corresponding fitness value of xopt
%
%   Author: Koen van der Blom, Hao Wang
%   Last modified: February 3, 2016
    
    % ----------------- general setting of variables ----------------------
    
    % TODO
    % static variables
    if dim == 2
        n = 12;   % the size of the problem: 12 for square and 9 for cube
    elseif dim == 3
        n = 9;    % the size of the problem: 12 for square and 9 for cube
    end
    
    % check input arguments
    if (dim ~= 2 && dim ~= 3)
        error('Invalid number of dimensions, use 2 or 3;');
    end
    
    % the pheno type of the solution is the permutation of integers 
    % from 1 to n ^ dim
    pheno_len = int16(n ^ dim);
    
    % TODO
    % At this point, you should think which geno type representation you
    % would like to use and thus determine the length of the geno type 
    % solution vector. Example of geno encoding: bit-string converted from 
    % array of integers. This is up to you :)
    geno_len = pheno_len;
    
    
    mu = 5;
    
    % TODO
    % endogenous parameters setting
    pc = 0.5;               % crossover rate
    pm = 0.8;               % mutation rate
    
    % plotting toggle
    is_plot = true;

	% internal counter variable
	evalcount = 0;     % count function evaluations
	gencount = 0;      % count generation/iterations
    
    % historical information
    hist_fitness = zeros(1, eval_budget);
    
    % ------------------- population initialization -----------------------
    % column vector representation is used throughout this code
    % you need to keep pheno type population updated with the geno types
    % for function evaluation
    
    % population
    pop_pheno = zeros(mu,pheno_len);            % pheno type
    pop_geno = zeros(mu,geno_len);              % geno type     
    fitness = zeros(mu, 1);                      % fitness values
    
    % TODO
    for i = 1 : mu
        pop_pheno(i , :) = randperm(pheno_len,pheno_len) ; 
 % generate pheno type individual uniformly
        pop_geno(i, :) = pop_pheno(i,:);  % convert them to geno type solution
        fitness(i) = fitness_func(pop_pheno(i,:));  % and evaluate the 
                                                     % solution...   
    end
    
    [fopt, index] = min(fitness);
    xopt = pop_geno(index,:);
    xopt_pheno = pop_pheno(index,:);
   
    % increase the evalcount by mu
    hist_fitness(evalcount+1:evalcount+mu) = fopt;
    evalcount = evalcount + mu;
    lambda = 5 * mu;
	% ----------------------- Evolution loop ------------------------------
	while evalcount < eval_budget
        % generate the a new population using crossover and mutation
%         pop_new_geno = zeros(mu,geno_len);
		for i = 1 : lambda % this mu here means the number of the children is the same with the parents

			% TODO
            % implement the selection operator.
            p_selection = zeros(mu,1);
            for iii  = 1 : mu 
                p_selection(iii) = fitness(iii) ./ sum(fitness);
            end
%             p_selection = fitness ./ sum(fitness);
%             size(p_selection')
%             size(pop_geno)
            rank_pop_pheno = [pop_geno,p_selection];
            
            j = 1;
            %select parent according to the fitness
            for ii = 1 : mu
                 if rank_pop_pheno(ii,geno_len+1) > (rand()/mu)
                        p(j,:) = pop_geno(ii,:);             % select parent from pop_geno according to fitness
                        j = j + 1;
                 end
            end
            if j > 1
				% TODO
                % implement the crossover operator
                selection = randperm(j-1,j-1);
                selection_index = selection(1);
            else
                selection_index = 1;
                p = pop_geno(1,:);
            end
            
			if (rand() < pc)
                pop_new_geno (i,:) = selfstrMutation(p(selection_index,:),geno_len); %self crossover
                for kk = 1 : 3
                    pop_new_geno(i,:) = selfp2pMutation(pop_new_geno(i,:),geno_len) ;    %because such crossover does not change the sequence of the candidatethus we use mutation to make it more random
                end
            else
                % No crossover, copy the parent chromosome
                pop_new_geno(i,:) =  p(selection_index,:);
                
            end
            
            % TODO
            % implement the muation operator
            if (rand() < pm)
            pop_new_geno(i,:) = selfp2pMutation(pop_new_geno(i,:),geno_len) ;      % apply the mutation and then  
                                             % store it in pop_new_geno
            end
            % TODO
            % repair the newly generated solution (if you want...)
            % the solution might be invalid because of duplicated integers
        end
        
        % Replace old population by the newly generated population
%         pop_geno = pop_new_geno;
       
        % TODO
        fitness_new = zeros(lambda,1);
        for i = 1 : lambda
            fitness_new(i) = fitness_func(pop_new_geno(i,:));     
        end
        pop_whole = [pop_geno;pop_new_geno];
        fitness_whole = [fitness;fitness_new];
        pop_whole_rank = [pop_whole,fitness_whole];
        pop_whole_rank = sortrows(pop_whole_rank,geno_len+1);
        
        pop_geno = pop_whole_rank(1:mu,1:geno_len);
        fitness = pop_whole_rank(1:mu,geno_len+1);
        % optimal solution in each iteration
%         [fopt_curr_iter, index] = min(fitness);
        fopt_curr_iter = fitness(1); 
        x_opt_curr_iter = pop_geno(1,:);
        x_opt_pheno_curr_iter = pop_geno(1,:);
        
        % keep track of the best solution ever found
        if fopt_curr_iter < fopt
            fopt = fopt_curr_iter;
            xopt = x_opt_curr_iter;
            xopt_pheno = x_opt_pheno_curr_iter;
        end
        
        % record historical information
        hist_fitness(evalcount+1:evalcount+mu) = fopt_curr_iter;
        
        % internal counters increment
		gencount = gencount + 1;
        evalcount = evalcount + mu;
        
        % Plot statistics
    end  
    
    
        if is_plot
            clf;
            subplot(1, 2, 1);
            plot(hist_fitness(1:evalcount));
            title('minimal error in the current iteration')
            
            subplot(1, 2, 2);
            bar(1:geno_len, xopt);
            xlim([1 geno_len]);
            title('best chromosome in the current iteration')
%             drawnow();
        end


end



%--------------------permutation for 2 dimension -------------------------------

% function s=permutation2d4(s,n,pm)
% %author: Li Yicheng
% %time :2016/3/13
% %s=permutation2d3(s,n,pm)   version3
% %input: s : candidate vector
% %       n : side of the square
% %       pm : the probability of mutation instead of croosover
% %output: modified vector
% 
% 
% s= reshape(s,[n,n]);
% std=  n * (n ^ 2 + 1) / 2;
% % ------------------------------------------------------------------
% %row
% P_row = ones(1,n);
% for i=1:n
%     if sum(s(i,:))==std
%         P_row(i)=0.125;
%     else if (abs(sum(s(i,:))-std)) < std*0.005
%             P_row(i)=0.80;
%         else if abs(sum(s(i,:))-std) <std*0.01
%              P_row(i)=0.875;
%             else if abs(sum(s(i,:))-std) <std*0.02
%                     P_row(i)=0.95;
%                 end
%             end
%         end
%     end
% end
% 
% 
% k=0;
% for i=1:n
%     if P_row(i) > rand(1)
%         index(k+1) = i;
%         k=k+1;
%     end
% end
% 
% % i_length = length(index);
% % index = index(randperm(i_length,i_length));
% 
% if (k >= 4)
%     index = index(randperm(k,k));
%     if pm>rand(1)
%          [s(index(1),:),s(index(2),:)]=p2p_mutation(s(index(1),:),s(index(2),:),n);
%          [s(:,index(3)),s(:,index(4))]=p2p_mutation(s(:,index(3)),s(:,index(4)),n);
%     else
%          [s(index(1),:),s(index(2),:)]=str_mutation(s(index(1),:),s(index(2),:),n);
%          [s(:,index(3)),s(:,index(4))]=str_mutation(s(:,index(3)),s(:,index(4)),n);
%     end
% else if (k>=2)
%     index = index(randperm(k,k));
%     if pm>rand(1)
%         [s(:,index(1)),s(:,index(2))]=p2p_mutation(s(:,index(1)),s(:,index(2)),n);
%     else
%         [s(:,index(1)),s(:,index(2))]=str_mutation(s(:,index(1)),s(:,index(2)),n);
%     end
%     else
%     disp('attention')
%     end
% end
% % ------------------------------------------------------------------
% %col
% % P_col = ones(1,n);
% % for i=1:n
% %     if sum(s(:,i))==std
% %         P_row(i)=0.125;
% %     else if (sum(s(:,i))-std) <std*0.005
% %             P_row(i)=0.80;
% %         else if (sum(s(:,i))-std) <std*0.01
% %              P_row(i)=0.875;
% %             else if (sum(s(:,i))-std) <std*0.02
% %                     P_row(i)=0.95;
% %                 end
% %             end
% %         end
% %     end
% % end
% % k=0;
% % for i=1:n
% %     if P_col(i) > rand(1)
% %         index(k+1) = i;
% %         k=k+1;
% %     end
% % end
% 
% % i_length = length(index);
% % index = index(randperm(i_length,i_length));
% % 
% % if (k >= 2)
% %     index = index(randperm(k,k));
% % 
% %     if pm>rand(1)
% %         [s(:,index(1)),s(:,index(2))]=p2p_mutation(s(:,index(1)),s(:,index(2)),n);
% %     else
% %         [s(:,index(1)),s(:,index(2))]=str_mutation(s(:,index(1)),s(:,index(2)),n);
% %     end
% % end
% 
% % ----------------------------------------------------------------------
% s=reshape(s,[1,n^2]);
% end
% 
% %-----------------------------------------------------------------------
% 
% 
% %--------------------permutation for 3 dimension-----------------------------
% function s=permutation3d6(s,n,pm)
% %author : Li Yicheng
% %date: 2016/3/13
% %s = permutation3d4(s,n,pm)   --version4
% % input : 
% %     s : candidate vector
% %     n : order of the problem
% %     pm : the probability of p2p mutation
% % output:
% %     a new candidate
% s_ori = s;
% s= reshape(s,[n,n,n]);
% std = n*(n^3 +1)/2; 
% 
% %-----------------------------------------------------------------
% % row
% 
% P_row=ones(1,n^2);
% for i = 1:n
%     for j=1:n
%         if abs(sum(s(j,:,i)))==std
%             P_row(j+(i-1)*n)=0.125;
%         else if abs(sum(s(j,:,i))-std)< std*0.001
%                   P_row(j+(i-1)*n)=0.5;
%             else if abs(sum(s(j,:,i))-std) < std*0.005
%                      P_row(j+(i-1)*n)=0.7;
%                 else if abs(sum(s(j,:,i))-std) < std*0.01
%                          P_row(j+(i-1)*n)=0.85;
%                     else if abs(sum(s(j,:,i))-std) < std*0.015
%                              P_row(j+(i-1)*n)=0.95;
%                         end
%                     end
%                 end
%             end
%         end
%      end
% end
% 
% k=0;
% for i=1:n^2
%     if P_row(i) > rand(1)
%         index(k+1) = i;
%         k= k+1;
%     end
% end
% if k>=4 %pick four as the parameter to mutate
%     index = index(randperm(k,k));
%     if mod(index(1),n)==0
%         j1=n;
%     else 
%         j1=mod(index(1),n);
%     end
%     if mod(index(2),n)==0
%         j2=n;
%     else 
%         j2=mod(index(2),n);
%     end
%     if mod(index(3),n)==0
%         j3=n;
%     else 
%         j3=mod(index(3),n);
%     end
%     if mod(index(4),n)==0
%         j4=n;
%     else 
%         j4=mod(index(4),n);
%     end
%     if pm>rand(1)
%          [s(ceil(index(1)/n),j1,:),...
%             s(ceil(index(2)/n),j2,:)]=...
%             p2p_mutation(...
%                 s(ceil(index(1)/n),j1,:),...
%                 s(ceil(index(2)/n),j2,:),...
%                 n);  
%          [s(:,j1,ceil(index(1)/n)),...
%             s(:,j2,ceil(index(2)/n))]=...
%             p2p_mutation(...
%                 s(:,j1,ceil(index(1)/n)),...
%                 s(:,j2,ceil(index(2)/n)),...
%                 n);
%          [s(j1,:,ceil(index(1)/n)),...
%             s(j2,:,ceil(index(2)/n))]=...
%             p2p_mutation(...
%                 s(j1,:,ceil(index(1)/n)),...
%                 s(j2,:,ceil(index(2)/n)),...
%                 n);
%          [s(ceil(index(3)/n),j3,:),...
%             s(ceil(index(4)/n),j4,:)]=...
%             p2p_mutation(...
%                 s(ceil(index(3)/n),j3,:),...
%                 s(ceil(index(4)/n),j4,:),...
%                 n);  
%          [s(:,j3,ceil(index(3)/n)),...
%             s(:,j4,ceil(index(4)/n))]=...
%             p2p_mutation(...
%                 s(:,j3,ceil(index(3)/n)),...
%                 s(:,j4,ceil(index(4)/n)),...
%                 n);
%          [s(j3,:,ceil(index(3)/n)),...
%             s(j4,:,ceil(index(4)/n))]=...
%             p2p_mutation(...
%                 s(j3,:,ceil(index(3)/n)),...
%                 s(j4,:,ceil(index(4)/n)),...
%                 n);
%     else
%          [s(ceil(index(1)/n),j2,:),...
%             s(ceil(index(2)/n),j1,:)]=...
%             str_mutation(...
%                 s(ceil(index(1)/n),j2,:),...
%                 s(ceil(index(2)/n),j1,:),...
%                 n);
%          [s(:,j1,ceil(index(1)/n)),...
%             s(:,j2,ceil(index(2)/n))]=...
%             str_mutation(...
%                 s(:,j1,ceil(index(1)/n)),...
%                 s(:,j2,ceil(index(2)/n)),...
%                 n);
%          [s(j1,:,ceil(index(1)/n)),...
%             s(j2,:,ceil(index(2)/n))]=...
%             str_mutation(...
%                 s(j1,:,ceil(index(1)/n)),...
%                 s(j2,:,ceil(index(2)/n)),...
%                 n);
%          [s(ceil(index(3)/n),j3,:),...
%             s(ceil(index(4)/n),j4,:)]=...
%             str_mutation(...
%                 s(ceil(index(3)/n),j3,:),...
%                 s(ceil(index(4)/n),j4,:),...
%                 n);  
%          [s(:,j3,ceil(index(3)/n)),...
%             s(:,j4,ceil(index(4)/n))]=...
%             str_mutation(...
%                 s(:,j3,ceil(index(3)/n)),...
%                 s(:,j4,ceil(index(4)/n)),...
%                 n);
%          [s(j3,:,ceil(index(3)/n)),...
%             s(j4,:,ceil(index(4)/n))]=...
%             str_mutation(...
%                 s(j3,:,ceil(index(3)/n)),...
%                 s(j4,:,ceil(index(4)/n)),...
%                 n);
%     end
%     s = reshape(s,[1,n^3]);
% else
%     disp('attention')
%     a = randperm(n^3,n^3);
%     [s_ori(a(1)),s_ori(a(2))] = p2p_mutation(s_ori(a(1)),s_ori(a(2)),n);
%     s = s_ori;
% end
% 
% end
% %-------------------------------------------------------------------------
% 
% %-------------assistance function-----------------------------------------
% function [a,b]=p2p_mutation(a,b,n)
% if ( a==b )
%     a = selfp2pMutation(a,n);
%     b = a;
% else
% mutate_index1 = round(rand(1)*(n-1)) +1;
% mutate_index2 = round(rand(1)*(n-1)) +1;
% a_mutation = a(mutate_index1);
% b_mutation = b(mutate_index2);
% a(mutate_index1)=b_mutation;
% b(mutate_index2)=a_mutation;
% end
% 
% end
% 
% function [a,b] = str_mutation(a,b,n)
% %input two vector each with n elements
% 
% length = round(rand(1)*(n-1))+1;
% index1=randperm(n+1-length);
% index2=randperm(n+1-length);
% a_part = a(index1(1):index1(1)+length -1);
% b_part = b(index2(1):index2(1)+length -1);
% a(index1:index1(1)+length - 1)=b_part;
% b(index2:index2(1)+length - 1 )=a_part;
% a = selfp2pMutation(a,n);
% b = selfp2pMutation(b,n);
% 
% end
% 
function  a = selfstrMutation(a,n)
length = round(rand(1)*(n-1))+1;
index1=int16(randperm(n+1-length,n+1-length));
index2=int16(randperm(n+1-length,n+1-length));
a_part = a(index1(1):index1(1)+length -1);
b_part = a(index2(1):index2(1)+length -1);
a(index1(1):(index1(1)+length - 1))=b_part;
b(index2(1):(index2(1)+length - 1 ))=a_part;
a = selfp2pMutation(a,n);
end
% 
% 
function a = selfp2pMutation(a,n)

mutate_index1 = round(rand(1)*(n-1)) +1;
mutate_index2 = round(rand(1)*(n-1)) +1;

a_mutation1 = a(mutate_index1);
a_mutation2 = a(mutate_index2);
a(mutate_index1)=a_mutation2;
a(mutate_index2)=a_mutation1;
end
