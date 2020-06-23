function [ E ] = calcEnergy( S2, S2_target)


    
    E = sum( (S2 - S2_target).^2 );
    
    

end

