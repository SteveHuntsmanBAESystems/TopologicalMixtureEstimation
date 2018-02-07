function p = unimixmaxjsd(x,p)

% Successively perturb a (presumed) unimodal mixture p(x) to obtain a new
% unimodal mixture with the same number of components but maximum
% Jensen-Shannon divergence.
%
% Copyright (c) 2018, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

%% Preliminaries
[M,n] = size(p);
JS = zeros(1,M*n);
JS(1) = jsd(x,p);

%% This loop will terminate if as assumed p(x) is a unimodal mixture
for k = 2:numel(JS)
    dp = zeros(M,M,n);
    obj = -Inf(M,M,n);
    for i = 1:M	% give to i
        ind_i = find(p(i,:)==max(p(i,:)));
        for j = 1:M	% take from j
            for r = 4:(n-3)
                %% Find saturating unimodality-preserving perturbation dp
                take = p(j,r)-min(p(j,r-1),p(j,r+1));
                if ismember(r,[ind_i(1)-1,ind_i,ind_i(end)+1])
                    give = Inf;
                else
                    give = max(p(i,r-1),p(i,r+1))-p(i,r);
                end
                dp(i,j,r) = min(take,give); 
                %% Evaluate jsd of a local test perturbation
                test = p;
                test(i,r) = test(i,r)+dp(i,j,r);
                test(j,r) = test(j,r)-dp(i,j,r);
                obj(i,j,r) = jsd(x,test);
            end
        end
    end
    %% Effect the test perturbation that leads to the largest jsd increase
    [JS(k),ind_max] = max(obj(:));
    [I,J,R] = ind2sub([M,M,n],ind_max);
    p(I,R) = p(I,R)+dp(I,J,R);
    p(J,R) = p(J,R)-dp(I,J,R);
    %% Terminate loop when further perturbation has no effect
    if dp(I,J,R) < max(p(:))*sqrt(eps), break; end % else dp(I,J,R); end
    if JS(k) <= JS(k-1)*(1+sqrt(eps)), break; end
end