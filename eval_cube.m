function f = eval_cube(x)
%   fitness function of the perfect magic cube: this code takes into
%   account the diagonal sums on each square slice
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
%           the mean squared error (MSE) of all each row, column, diagonal
%           and space diagonal sum to the magic constant is computed
%
%   Author: Koen van der Blom, Hao Wang
%   Last modified: February 3, 2016

%     n = length(x) ^ (1 / 3);
    n = nthroot(length(x),3);
    
    if round(n) ~= n
        error(['Invalid length of the solution! The length of solutions ',...
            'should be cubic numbers']);
    end
    
    magic_constant = n * (n ^ 3 + 1) / 2;
    
    % convert the vector representation into a cube (3-D matrix)
    X = reshape(x, [n, n, n]);
    
    f = 0;
    counter = 0;
    
    % compute all the column sums and row sums
    for i = 1 : n
        
        square = X(:, :, i);
       
        col_sum = sum(square, 1);
        row_sum = sum(square, 2);
        diag_sum = [sum(diag(square)), sum(diag(rot90(square)))];
        counter = counter + 2 * n + 2;
%     
%         [col_sum, row_sum', diag_sum]
%         sum(([col_sum, row_sum', diag_sum] - magic_constant) .^ 2)
        f = f + sum(([col_sum, row_sum', diag_sum] - magic_constant) .^ 2);
        
    end

    % compute all the pillar sum
    for i = 1 : n
        square = reshape(X(:, i, :), [n, n]);
        pillar_sum = sum(square, 2);
        diag_sum = [sum(diag(square)), sum(diag(rot90(square)))];
        
        counter = counter + n + 2;

        f = f + sum(([pillar_sum', diag_sum] - magic_constant) .^ 2);
    end
    
    % the rest diagonals
    for i = 1 : n
        square = reshape(X(i, :, :), [n, n]);
        diag_sum = [sum(diag(square)), sum(diag(rot90(square)))];
        
        counter = counter + 2;

        f = f + sum((diag_sum - magic_constant) .^ 2);
    end
    
    % four space diagonals
    space_square1 = zeros(n,n);
    for i=1:n
        space_square1(i,:)=diag(X(:,:,i))';
    end
    space_square2 = zeros(n,n);
    for i=1:n
        space_square2(i,:)=diag(rot90(X(:,:,i)))';
    end
   
    
    space_diag_sum = [sum(diag(space_square1)), sum(diag(rot90(space_square1))),...
        sum(diag(space_square2)), sum(diag(rot90(space_square2)))];
    
    f = f + sum((space_diag_sum - magic_constant) .^ 2);
    
    counter = counter + 4;
    
    f = f / counter;
    
end