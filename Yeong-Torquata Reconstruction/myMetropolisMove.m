function [move, p, R] = myMetropolisMove(dE, T)
    
    % Generate random number between 0 and 1
    R = rand(1);
    % Calculate acceptance probability
    p = exp(-dE/T);
    % if the energy change is less than 0, accept the move automatically
    if dE <= 0
        move = 1;
    % if dE > 0, accept with probability exp(-dE/T)
    elseif R <= p
        move = 1;
    else 
        move = 0;
    end
    
end