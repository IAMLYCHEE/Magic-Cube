function f = eval_semi_cube(x)
%   fitness function of the semi-perfect magic cube: this code does not 
%   takes into account the diagonal sums on each square slice
%
%   Parameters
%   ----------
%       x : array, the solution vector that represents a magic cube.
%           By default, the solution vector is converted to a magic cube
%           vertically from the top square to the bottom square and each
%           square is constructed columnwisely.
%           
%   Output
%   ----------
%       f : double, the error value of the input solution vector.
%           the mean squared error (MSE) of all each row, column and
%           space diagonal sum to the magic constant is computed
%
%   Author: Koen van der Blom, Hao Wang
%   Last modified: February 3, 2016


    n = round(length(x) ^ (1 / 3));
    
%     if round(n) ~= n
%         error('Invalid length of the solution!');
%     end
    
    magic_constant = n * (n ^ 3 + 1) / 2;
    
    % the cube...
    X = reshape(x, [n, n, n]);
    
    f = 0;
    counter = 0;
    
    % compute all the column sums and row sums
    for i = 1 : n
        
        square = X(:, :, i);
       
        col_sum = sum(square, 1);
        row_sum = sum(square, 2);
        counter = counter + 2 * n;
    
        f = f + sum(([col_sum, row_sum'] - magic_constant) .^ 2);
        
    end

    % compute all the pillar sum
    for i = 1 : n
        square = reshape(X(:, i, :), [n, n]);
        pillar_sum = sum(square, 2);
        
        counter = counter + n;

        f = f + sum((pillar_sum - magic_constant) .^ 2);
    end
    
    % four space diagonals
    space_square1 = [diag(X(:, :, 1)), diag(X(:, :, 2)), diag(X(:, :, 3))];
    space_square2 = [diag(rot90(X(:, :, 1))), diag(rot90(X(:, :, 2))), diag(rot90(X(:, :, 3)))];
    
    space_diag_sum = [sum(diag(space_square1)), sum(diag(rot90(space_square1))),...
        sum(diag(space_square2)), sum(diag(rot90(space_square2)))];
    
    f = f + sum((space_diag_sum - magic_constant) .^ 2);
    
    counter = counter + 4;
    
    f = f / counter;
    
    % Update run statistics
    statistics(f);
    
end