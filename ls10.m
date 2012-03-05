function x = ls10(low, high, steps)

% x = ls10(low, high, steps)
%
% Shortcut for logspace(log10(low), log10(high), steps)

x = logspace(log10(low), log10(high), steps);
