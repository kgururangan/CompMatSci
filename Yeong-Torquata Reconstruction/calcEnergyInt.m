function [ E ] = calcEnergyInt( S2, S2_target,rsamp)

    
    % Create domain vector for integration in polar coordinates
    r = 0:rsamp;
    r = r';
   
    % Define the integrand
    I = 2*pi*r.*(S2-S2_target).^2;
    
    % E is the integral over polar coordinates of the squared difference between correlation
    % function
    E = trapz(r, I);
    

    
    

end

