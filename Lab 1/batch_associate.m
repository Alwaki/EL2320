% This function performs the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           z(t)                2Xn
% Outputs: 
%           c(t)                1Xn
%           outlier             1Xn
%           nu_bar(t)           2nX1
%           H_bar(t)            2nX3
function [c, outlier, nu_bar, H_bar] = batch_associate(mu_bar, sigma_bar, z)
        
        % YOUR IMPLEMENTATION %
        n = size(z, 2);
        H_bar = zeros(1, 6*n);
        nu_bar = zeros(2, n);
        outlier = zeros(1, n);
        c = zeros(1, n);
        
        for i = 1:n
            j = (i*6) - 5;
            
            z_i = z(:,i);
            [c_i,outlier_i, nu, S, H] = associate(mu_bar, sigma_bar, z_i);          
            H_i = H(:,:,c_i);
            nu_i = nu(:,c_i);
            
            c(:,i) = c_i;            
            outlier(:,i) = outlier_i;
            nu_bar(:,i) = nu_i;
            H_bar(j:i*6) = reshape(H_i', 1, 6);
   
        end
        
        nu_bar = nu_bar(:);
        H_bar = reshape(H_bar, 3, 2*n)';

end