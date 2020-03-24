classdef mapsi_data
%MAPSI_DATA creates a data object used to pass all of the collected data
% and relevant instrument parameters to the algorithm.
%
% Usage: data = mapsi_data(Q, I, sigma, c, b, kernel)
%
% INPUTS:
% Q - wave vectors at which intensities are  collected (an Nw x 3 matrix)
% I - measured intensities
% sigma - measured errors
% c - material proportionality constant
% b - baseline
% kernel - function of two arguments (in and out) to act as the kernel for
%   the Fredholm integral
%
% OUTPUTS:
% data - a mapsi_data object

properties
    Q
    I
    sigma
    c
    b
    kernel
end

methods
    function data = mapsi_data( Q, I, sigma, c, b, kernel)
        if nargin ~= 6
            error('An incorrect number of arguments was provided')
        end
        
        data.Q = Q;
        data.I = I;
        data.sigma = sigma;
        data.c =c;
        data.b = b;
        data.kernel = kernel;
    end
end
    
end