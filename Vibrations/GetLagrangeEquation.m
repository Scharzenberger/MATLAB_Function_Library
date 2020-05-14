function Leqs = GetLagrangeEquation(KE, DP, PE, vars, dvars, ddvars)

% This function computes the Lagrange equation for a mass-spring-damper system.

% Retrieve the number of degrees of freedom.
DOF = length(vars);

% Preallocate an array to store the Lagrange equations.
Leqs = sym(zeros(DOF, 1));

% Compute each Lagrange equation.
for k = 1:DOF       % Iterate through each DOF...
    
    Leqs(k) = subs(diff(KE, dvars(k)), dvars(k), ddvars(k)) - diff(KE, vars(k)) + diff(PE, vars(k)) + diff(DP, dvars(k))
    
    
end




end

