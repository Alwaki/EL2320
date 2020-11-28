% This function performs multinomial re-sampling
% Inputs:   
%           S_bar(t):       4XM
% Outputs:
%           S(t):           4XM
function S = multinomial_resample(S_bar)

    global M % number of particles
    
    % YOUR IMPLEMENTATION
    CDF = cumsum(S_bar(4, :));
    S = zeros(4, M);
    
    for m = 1:M
        rm = rand;
        xi = find(CDF >= rm, 1);
        S(1:3, m) = S_bar(1:3, xi);
    end
    
    S(4,:) = 1/M * ones(size(S_bar(4,:)));
        
end
