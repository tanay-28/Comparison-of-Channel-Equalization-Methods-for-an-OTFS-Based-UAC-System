function [cdf, pdf_values] = pdf_channel(r,K,mu,sigma)
%this function computes the cdf from the given values 
% Define the range of 'h' values
h_max = 1; % Adjust the maximum 'h' value as needed
num_points = 1000;
h = linspace(0.01, h_max, num_points);
% h = h_max/num_points*(1:num_points);
% h = 0.001*(1:num_points);
%% Initialize an array to store PDF values
pdf_values = zeros(size(h));
cdf=zeros(size(h));
% Loop over 'h' values
for i = 1:length(h)
 
    sum=0;
    for S=0.01:0.01:10
    h_i = h(i);
    p = (2 * (K + 1) * h_i / (sqrt(2 * pi) * r * sigma * S^3)) ...
        * exp(-((K + 1) * h_i^2 / S*S + K + ((log(S) - mu) / (sqrt(2) * r * sigma))^2)) ...
        * besselk(0, 2 * h_i / sqrt(K * (K + 1)));
    sum=sum+p;
    end
    %% Numerical integration using the quad function
    pdf_values(i) = sum; % Adjust the upper limit as needed
    if(i>1)
    cdf(i)=cdf(i-1)+pdf_values(i);
    end
end    
% pdf_rms=rms(pdf_values);
pdf_values=pdf_values/norm(pdf_values);
% Normalize the CDF to make it monotonically increasing and between 0 and 1
cdf = cdf / cdf(end);
end