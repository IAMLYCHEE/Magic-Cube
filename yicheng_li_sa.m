function [xopt, fopt,hist_best_f] = yicheng_li_sa(dim, eval_budget, fitness_func)
% [xopt, fopt] = sa_skeleton(dim, eval_budget, fitness_func)
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

	% Set true to do online statistics plotting
	is_plot = false;
    
    % TODO
	% Initialize static parameters
% 	pm = ...         % mutation rate
% 	alpha = ...      % temperature decaying parameter
% 	k = ...
%     
% I use the logarithmic as cooling schedule
    if dim == 2 
        k = 5;
    else 
        k = 15;
    end
    pm_init=0.1;
    pm=pm_init;
    
    if dim == 2
        len = 12^2;    % length of the solution vector, shoud be 12^2 
                     % when dim == 2
    elseif dim == 3
        len = 9^3;    % length of the solution vector, shoud be 9^3 when 
                     % dim == 3
    end
    
    % TODO
    % Set initial temperature
    if(dim == 3)
    	T_init = 6000;
    
    else
        T_init=500;
    end
    
    T_para=2;
    T=T_init;

	% Statistics data
	evalcount = 0;
	gencount = 0;
	fopt = Inf;
	xopt = NaN(len, 1);
	hist_best_f = NaN(1, eval_budget);
	hist_f = NaN(1, ceil(eval_budget / k));
	hist_temperature = NaN(1, ceil(eval_budget / k));
     
    % TODO
	% Generate initial solution and evaluate
	s = randperm(len,len) ; 
%     s = x_init;
	f = fitness_func(s);         % evaluate the solution using fitness_func
	
	if (f < fopt)
		fopt = f;
		xopt = s;
    end
    

    evalcount = evalcount + 1;             % Increase evaluation counter
	hist_best_f(evalcount) = fopt;       
    p2pCounter = 0;
    strCounter = 0;
    % evolution loop
    
    local_count = 0;
	while evalcount < eval_budget

		for j = 1 : k
            
%             % TODO
% 			s_new = ...   % Generate a new solution by the permutation of s
% 			f_new = ...   % evaluate the new solution
            if(dim == 3)
                s_new = permutation3d6(s,nthroot(len,3),pm);
                f_new = fitness_func(s_new);
            else
                s_new = permutation2d4(s,sqrt(len),pm);
                f_new = fitness_func(s_new);
            end
               
            
			evalcount = evalcount + 1;   % Increase evaluation counter 
			gencount = gencount + 1;     % Increase iteration counter
			
			if (f_new < f)
                % TODO
                % accept the new solution when it is better
                s=s_new;
                f=f_new;
%                 local_count = 0;
            else
                % TODO
                % choose to accept or reject the new solution 
                % probabilistically based on the current temperature
                p=exp(-(f_new-f)/(T));
                if p<0.0001 
                    p=0.0001; 
                end
                if(p>rand(1))
                    s=s_new;
                    f=f_new;
                end
%                 local_count = local_count + 1
             end
        
            % update the best solution found so far
			if (f < fopt)
				fopt = f;
				xopt = s;
%                 local_count = 0;
%             else 
%                 local_count = local_count + 1
            end
            
            if (fopt<50000)
                pc=0.5;
            end
            
            hist_best_f(evalcount) = fopt;   % tracking the best fitness
                                             % ever found

			% Generation best statistics
			hist_f(gencount) = f;
			hist_temperature(gencount) = T;

			if (is_plot)
				% Plot statistics
				clf

				subplot(1, 3, 1)
				plot(hist_best_f(evalcount-2000:evalcount))
				title('minimal error found so far')

				subplot(1, 3, 2)
				plot(hist_f(evalcount - 2000:gencount))
				title('minimal error in the current generation')

				subplot(1, 3, 3)
				plot(hist_temperature(1:gencount))
				title('Temperature')

				drawnow()

			end

        end
%         fopt;
        % TODO
		% Temperature update
%         p
        T_para= T_para +0.05;
		T = T_init/T_para;
        pm = 1/log(exp(1)+T_para) + 0.6;
        
%         if(local_count >= 200)
%             T = T_init;
%             T_para = 2;
%             local_count = 0;
%         end
%         if pm > rand(1)
%             p2pCounter = p2pCounter + 1;
%         else
%             strCounter = strCounter + 1;
%         end
%         pm=1;
    end
    
    
    			subplot(1, 3, 1)
				plot(hist_best_f)
				title('minimal error found so far')

				subplot(1, 3, 2)
				plot(hist_f)
				title('minimal error in the current generation')

				subplot(1, 3, 3)
				plot(hist_temperature(1:gencount))
				title('Temperature')
                
%                 p2pCounter;
%                 strCounter;
%                 operation = p2pCounter / (strCounter + p2pCounter)

end




%--------------------permutation for 2 dimension -------------------------------

function s=permutation2d4(s,n,pm)
%author: Li Yicheng
%time :2016/3/13
%s=permutation2d3(s,n,pm)   version3
%input: s : candidate vector
%       n : side of the square
%       pm : the probability of mutation instead of croosover
%output: modified vector


s= reshape(s,[n,n]);
std=  n * (n ^ 2 + 1) / 2;
% ------------------------------------------------------------------
%row
P_row = ones(1,n);
for i=1:n
    if sum(s(i,:))==std
        P_row(i)=0.125;
    else if (abs(sum(s(i,:))-std)) < std*0.005
            P_row(i)=0.80;
        else if abs(sum(s(i,:))-std) <std*0.01
             P_row(i)=0.875;
            else if abs(sum(s(i,:))-std) <std*0.02
                    P_row(i)=0.95;
                end
            end
        end
    end
end


k=0;
for i=1:n
    if P_row(i) > rand(1)
        index(k+1) = i;
        k=k+1;
    end
end

% i_length = length(index);
% index = index(randperm(i_length,i_length));

if (k >= 4)
    index = index(randperm(k,k));
    if pm>rand(1)
         [s(index(1),:),s(index(2),:)]=p2p_mutation(s(index(1),:),s(index(2),:),n);
         [s(:,index(3)),s(:,index(4))]=p2p_mutation(s(:,index(3)),s(:,index(4)),n);
    else
         [s(index(1),:),s(index(2),:)]=str_mutation(s(index(1),:),s(index(2),:),n);
         [s(:,index(3)),s(:,index(4))]=str_mutation(s(:,index(3)),s(:,index(4)),n);
    end
else if (k>=2)
    index = index(randperm(k,k));
    if pm>rand(1)
        [s(:,index(1)),s(:,index(2))]=p2p_mutation(s(:,index(1)),s(:,index(2)),n);
    else
        [s(:,index(1)),s(:,index(2))]=str_mutation(s(:,index(1)),s(:,index(2)),n);
    end
    else
    disp('attention')
    end
end
% ------------------------------------------------------------------
%col
% P_col = ones(1,n);
% for i=1:n
%     if sum(s(:,i))==std
%         P_row(i)=0.125;
%     else if (sum(s(:,i))-std) <std*0.005
%             P_row(i)=0.80;
%         else if (sum(s(:,i))-std) <std*0.01
%              P_row(i)=0.875;
%             else if (sum(s(:,i))-std) <std*0.02
%                     P_row(i)=0.95;
%                 end
%             end
%         end
%     end
% end
% k=0;
% for i=1:n
%     if P_col(i) > rand(1)
%         index(k+1) = i;
%         k=k+1;
%     end
% end

% i_length = length(index);
% index = index(randperm(i_length,i_length));
% 
% if (k >= 2)
%     index = index(randperm(k,k));
% 
%     if pm>rand(1)
%         [s(:,index(1)),s(:,index(2))]=p2p_mutation(s(:,index(1)),s(:,index(2)),n);
%     else
%         [s(:,index(1)),s(:,index(2))]=str_mutation(s(:,index(1)),s(:,index(2)),n);
%     end
% end

% ----------------------------------------------------------------------
s=reshape(s,[1,n^2]);
end

%-----------------------------------------------------------------------


%--------------------permutation for 3 dimension-----------------------------
function s=permutation3d6(s,n,pm)
%author : Li Yicheng
%date: 2016/3/13
%s = permutation3d4(s,n,pm)   --version4
% input : 
%     s : candidate vector
%     n : order of the problem
%     pm : the probability of p2p mutation
% output:
%     a new candidate
s_ori = s;
s= reshape(s,[n,n,n]);
std = n*(n^3 +1)/2; 

%-----------------------------------------------------------------
% row

P_row=ones(1,n^2);
for i = 1:n
    for j=1:n
        if abs(sum(s(j,:,i)))==std
            P_row(j+(i-1)*n)=0.125;
        else if abs(sum(s(j,:,i))-std)< std*0.001
                  P_row(j+(i-1)*n)=0.5;
            else if abs(sum(s(j,:,i))-std) < std*0.005
                     P_row(j+(i-1)*n)=0.7;
                else if abs(sum(s(j,:,i))-std) < std*0.01
                         P_row(j+(i-1)*n)=0.85;
                    else if abs(sum(s(j,:,i))-std) < std*0.015
                             P_row(j+(i-1)*n)=0.95;
                        end
                    end
                end
            end
        end
     end
end

k=0;
for i=1:n^2
    if P_row(i) > rand(1)
        index(k+1) = i;
        k= k+1;
    end
end
if k>=4 %pick four as the parameter to mutate
    index = index(randperm(k,k));
    if mod(index(1),n)==0
        j1=n;
    else 
        j1=mod(index(1),n);
    end
    if mod(index(2),n)==0
        j2=n;
    else 
        j2=mod(index(2),n);
    end
    if mod(index(3),n)==0
        j3=n;
    else 
        j3=mod(index(3),n);
    end
    if mod(index(4),n)==0
        j4=n;
    else 
        j4=mod(index(4),n);
    end
    if pm>rand(1)
         [s(ceil(index(1)/n),j1,:),...
            s(ceil(index(2)/n),j2,:)]=...
            p2p_mutation(...
                s(ceil(index(1)/n),j1,:),...
                s(ceil(index(2)/n),j2,:),...
                n);  
         [s(:,j1,ceil(index(1)/n)),...
            s(:,j2,ceil(index(2)/n))]=...
            p2p_mutation(...
                s(:,j1,ceil(index(1)/n)),...
                s(:,j2,ceil(index(2)/n)),...
                n);
         [s(j1,:,ceil(index(1)/n)),...
            s(j2,:,ceil(index(2)/n))]=...
            p2p_mutation(...
                s(j1,:,ceil(index(1)/n)),...
                s(j2,:,ceil(index(2)/n)),...
                n);
         [s(ceil(index(3)/n),j3,:),...
            s(ceil(index(4)/n),j4,:)]=...
            p2p_mutation(...
                s(ceil(index(3)/n),j3,:),...
                s(ceil(index(4)/n),j4,:),...
                n);  
         [s(:,j3,ceil(index(3)/n)),...
            s(:,j4,ceil(index(4)/n))]=...
            p2p_mutation(...
                s(:,j3,ceil(index(3)/n)),...
                s(:,j4,ceil(index(4)/n)),...
                n);
         [s(j3,:,ceil(index(3)/n)),...
            s(j4,:,ceil(index(4)/n))]=...
            p2p_mutation(...
                s(j3,:,ceil(index(3)/n)),...
                s(j4,:,ceil(index(4)/n)),...
                n);
    else
         [s(ceil(index(1)/n),j2,:),...
            s(ceil(index(2)/n),j1,:)]=...
            str_mutation(...
                s(ceil(index(1)/n),j2,:),...
                s(ceil(index(2)/n),j1,:),...
                n);
         [s(:,j1,ceil(index(1)/n)),...
            s(:,j2,ceil(index(2)/n))]=...
            str_mutation(...
                s(:,j1,ceil(index(1)/n)),...
                s(:,j2,ceil(index(2)/n)),...
                n);
         [s(j1,:,ceil(index(1)/n)),...
            s(j2,:,ceil(index(2)/n))]=...
            str_mutation(...
                s(j1,:,ceil(index(1)/n)),...
                s(j2,:,ceil(index(2)/n)),...
                n);
         [s(ceil(index(3)/n),j3,:),...
            s(ceil(index(4)/n),j4,:)]=...
            str_mutation(...
                s(ceil(index(3)/n),j3,:),...
                s(ceil(index(4)/n),j4,:),...
                n);  
         [s(:,j3,ceil(index(3)/n)),...
            s(:,j4,ceil(index(4)/n))]=...
            str_mutation(...
                s(:,j3,ceil(index(3)/n)),...
                s(:,j4,ceil(index(4)/n)),...
                n);
         [s(j3,:,ceil(index(3)/n)),...
            s(j4,:,ceil(index(4)/n))]=...
            str_mutation(...
                s(j3,:,ceil(index(3)/n)),...
                s(j4,:,ceil(index(4)/n)),...
                n);
    end
    s = reshape(s,[1,n^3]);
else
    disp('attention')
    a = randperm(n^3,n^3);
    [s_ori(a(1)),s_ori(a(2))] = p2p_mutation(s_ori(a(1)),s_ori(a(2)),n);
    s = s_ori;
end

end
%-------------------------------------------------------------------------

%-------------assistance function-----------------------------------------
function [a,b]=p2p_mutation(a,b,n)
if ( a==b )
    a = selfp2pMutation(a,n);
    b = a;
else
mutate_index1 = round(rand(1)*(n-1)) +1;
mutate_index2 = round(rand(1)*(n-1)) +1;
a_mutation = a(mutate_index1);
b_mutation = b(mutate_index2);
a(mutate_index1)=b_mutation;
b(mutate_index2)=a_mutation;
end

end

function [a,b] = str_mutation(a,b,n)
%input two vector each with n elements

length = round(rand(1)*(n-1))+1;
index1=randperm(n+1-length);
index2=randperm(n+1-length);
a_part = a(index1(1):index1(1)+length -1);
b_part = b(index2(1):index2(1)+length -1);
a(index1:index1(1)+length - 1)=b_part;
b(index2:index2(1)+length - 1 )=a_part;
a = selfp2pMutation(a,n);
b = selfp2pMutation(b,n);

end

function a = selfp2pMutation(a,n)

mutate_index1 = round(rand(1)*(n-1)) +1;
mutate_index2 = round(rand(1)*(n-1)) +1;

a_mutation1 = a(mutate_index1);
a_mutation2 = a(mutate_index2);
a(mutate_index1)=a_mutation2;
a(mutate_index2)=a_mutation1;
end






