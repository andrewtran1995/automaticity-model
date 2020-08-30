function [idx] = rand_discrete(distr)
    % Given discrete distribution, return index of chosen index.
    cum_distr = cumsum(distr);
    idx = find(rand<cum_distr, 1);
end