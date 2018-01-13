function [ wmparm ] = autowmrnd( mu, range )
%AUTOWMRND Return randomly generated value for working memory parameter
%   Uses input arguments to randomly generate value for automaticity
%   model's working memory parameters. Random values distributed uniformly.
    wmparm = mu + (-range + 2*range*rand);
    wmparm = max(eps,wmparm);
end

