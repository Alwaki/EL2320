% This function performs the ML data association
%           S_bar(t)                 4XM
%           z(t)                     2Xn
%           association_ground_truth 1Xn | ground truth landmark ID for
%           every measurement  1 (actually not used)
% Outputs: 
%           outlier                  1Xn
%           Psi(t)                   1XnXM
%           c                        1xnxM (actually not ever used)
function [outlier, Psi, c] = associate(S_bar, z, association_ground_truth)
    if nargin < 3
        association_ground_truth = [];
    end

    global lambda_psi % threshold on average likelihood for outlier detection
    global Q % covariance matrix of the measurement model
    global M % number of particles
    global N % number of landmarks
    
    % YOUR IMPLEMENTATION
    n = size(z, 2);
    %Preallocate
    outlier = zeros(1, n);
    Psi = zeros(1, n, M);
    c = zeros(1, n, M);
    
    for i = 1:n
         for k = 1:N
             z_hat = observation_model(S_bar, k);
             nu = z(:, i) - z_hat;
             nu(2, :) = mod(nu(2, :) + pi, 2 * pi) - pi;
             psi = 1/(2 * pi * sqrt(det(Q))) * exp(-1/2 * nu' / Q * nu); %POSSIBLE ERROR IN SQRT
             psi = diag(psi);
             psi_k(:, k) = psi;
         end
        [Psi_max, ind] = max(psi_k, [], 2);
        c(1, i, :) = ind;
        Psi(1, i, :) = Psi_max;
        outlier(1, i) = (1/M) * sum(Psi_max) <= lambda_psi;
    end
end


