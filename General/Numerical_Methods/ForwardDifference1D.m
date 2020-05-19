function dfdx = ForwardDifference1D(fs, xs, n)

% This function computes the 1D nth order forward difference derivative approximation of fs given grid spacing xs.

% Inputs:
    % fs = Array of scalar function values.
    % xs = Array of scalar grid points.
    % n = Derivative Order.  i.e., n = 1 is the first derivative, n = 2 is the second derivative, etc.

% Outputs:   
    % dfdx = nth Derivative Approximation of fs given grid spacing xs.

% Retrieve the number of entries.
num_entries = length(xs);

% Preallocate a variable to store the forward difference approximations.
dfdx = zeros(num_entries - n, 1);

% Iterate through each of the entries.
for k = 1:(num_entries - n)             % Iterate through each of the entries for which we have enough information to compute a nth order forward difference...
    
    % Initialize the dx and df calculations.
    dx = 1;
    df = fs(k + n);
    
    % Compute the dx and df approximations via an nth order foward difference method.
    for i = 1:n             % Iterate through each of the orders...
        
        % Approximate the dx term.
        dx = dx.*(xs(k + i) - xs(k + i - 1));

        % Approximate the df term.
        df = df + ((-1).^i).*nchoosek(n, i).*fs(k + n - i);
        
    end
    
    % Approximate the nth order df/dx calculation.
    dfdx(k) = df./dx;
    
end

% Determine whether to transpose the output.
if size(xs, 2) > size(xs, 1)            % If the input has more columns than rows...
    
    % Transpose the output.
    dfdx = dfdx';
    
end

end

