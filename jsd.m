function y = jsd(x,p)

% (Differential) Jensen-Shannon divergence (in nats) of a continuous
% mixture p with argument x. p is ASSUMED to be appropriately normalized:
% i.e., the integral of sum(p,1) w.r.t x is assumed to be unity.
%
% Copyright (c) 2018, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

if size(p,2) ~= numel(x), error('must have size(p,2) == numel(x)'); end
pi = sum(p,2)/sum(p(:));
y = entropy(x,pi(:)'*p)-pi(:)'*entropy(x,p);

end

%% Local function

function S = entropy(x,p)

% Rows of p correspond to values of continuous PDFs at a common argument x.
% This computes the (differential) Shannon entropy of rows of p in nats.

dx = diff(x(:)');   % make dx a row array for extra safety w/ bsxfun
p = bsxfun(@rdivide,p,sum(bsxfun(@times,p,[mean(dx),dx]),2));   % normalize
S = -sum(bsxfun(@times,p.*log(p+(p==0)),[mean(dx),dx]),2);

end