function f = eval_square(x)
%   fitness function of the magic square
%
%   Parameters
%   ----------
%       x : array, the solution vector that represents a magic square.
%           By default, the solution vector is converted to a magic square
%           columnwisely.
%   Output
%   ----------
%       f : double, the error value of the input solution vector.
%           the mean squared error (MSE) of all each row, column and
%           diagonal sum to the magic constant is computed
%
%   Author: Koen van der Blom, Hao Wang
%   Last modified: February 3, 2016

    n = sqrt(length(x));
    
    if round(n) ~= n
         error(['Invalid length of the solution! The length of solutions ',...
            'should be square numbers']);
    end
    
    magic_constant = n * (n ^ 2 + 1) / 2;
    
     % convert the vector representation into a square (matrix)
    X = reshape(x, n, n);
   
    col_sum = sum(X, 1);
    row_sum = sum(X, 2)';
    diag_sum = [sum(diag(X)), sum(diag(rot90(X)))];
    
    f = mean(([col_sum, row_sum, diag_sum] - magic_constant) .^ 2);
    
end



