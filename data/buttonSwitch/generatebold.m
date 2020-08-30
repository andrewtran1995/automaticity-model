function [ BOLD ] = generatebold( stimvector )
%GENERATEBOLD Convolve stimulus vector w/ HRF for predicted BOLD
%   Convolve stimulus vector in increments of 0.01 seconds with the HRF
%   stimvector is expected to be in time increments of 0.001 seconds
    
    % Set parameter values of hrf
    t1 = 0.01; n = 4; lamda = 2;
    
    % Define time axis
    seconds = length(stimvector)/1000;
    t = 0.01:0.01:seconds;
    % Create hrf
    hrf = ((t-t1).^(n-1)).*exp(-(t-t1)/lamda)/((lamda^n)*factorial(n-1));
    
    % Compute convolution, dividing by 100 to set time unit at 0.01 seconds
    % Interpret stimvector in increments of 0.01 seconds
    BOLD = conv(stimvector(1:10:length(stimvector)), hrf')/100;
end

