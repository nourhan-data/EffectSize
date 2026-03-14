function [fit] = CEC2011func(x, fun, ~)
% Benchmark functions from the CEC 2011 competition
% Developed by Dan Simon, based on Mahamed Omran's CEC 2005 func.m
% INPUTS:
%   x is a D-element row vector representing the candidate solution we want to evalaute, where D is the problem dimension
%   fun is the index of the function of interest
%   opt_f is the minimum function value
% OUTPUT:
%   fit is the fitness of x (i.e. fit = f(x))
% NOTE: This files assumes that you include the following lines in your file:
%       global  initial_flag
%       initial_flag = 0;
global initial_flag
persistent MGADSMproblem
switch(fun)
    case {1, 2, 3, 4, 5, 6, 7, 8}
        fit = bench_func(x, fun);
    case {9}
        fit = cost_fn(x);
    case {10}
        fit = antennafunccircular(x, [50,120], 180, 0.5);
    case {111}
        fit = fn_DED_5(x);
    case {112}
        fit = fn_DED_10([x, 55*ones(1,24)]);
    case {113}
        fit = fn_ELD_6(x);
    case {114}
        fit = fn_ELD_13(x);
    case {115}
        fit = fn_ELD_15(x);
    case {116}
        fit = fn_ELD_40(x);
    case {117}
        fit = fn_ELD_140(x);
    case {118}
        fit = fn_HT_ELD_Case_1(x);
    case {119}
        fit = fn_HT_ELD_Case_2(x);
    case {1110}
        fit = fn_HT_ELD_Case_3(x);
    case {12}
        if initial_flag == 0
            load('messengerfull.mat');
            initialflag = 1;
        end
        fit = messengerfull(x, MGADSMproblem);
    case {13}
        if initial_flag == 0
            load('cassini2.mat');
            initialflag = 1;
        end
        fit = cassini2(x, MGADSMproblem);
    otherwise
        error('No such function');
end
end