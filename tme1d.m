function y = tme1d(X,plotflag)

% Topological mixture estimation using Gaussian kernel and
% deblurring/reblurring technique.
%
% If plotflag is set, the results are displayed;
% if plotflag > 0, a legend is also displayed on the plots
%
% Besides the outputs of tde1d, this also produces four mixtures, with
% respective fields sweep, unblurred, deblurred, and reblurred. Typically
% the last of these is the only desirable one, though all but the second
% have densities that agree with the results of tde1d.
%
% Copyright (c) 2018, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

%% Preliminaries
X = X(:);
n = numel(X);
nh = min(n,100);
x = linspace(min(X),max(X),nh);
DX = max(X)-min(X);
h = DX./(1:nh);

%% TDE
tde = tde1d(X,1);   % Gaussian kernel

%% Unblurred unimodal decomposition and JSD
sweep = unidec(tde.y,0);
M = size(sweep,1);
sweep = sweep(:,2:end-1);   % see plot shift below
unblurred = unimixmaxjsd(x,sweep);

%% Deblurred bandwidth and density estimate
h_minus = min(h(tde.uc==tde.mfuc));
delta_h = sqrt(tde.h^2-h_minus^2);
p_hat = zeros(1,numel(x));
for k = 1:n
    p_hat = p_hat+exp(-.5*((x-X(k))/h_minus).^2)/(h_minus*sqrt(2*pi));
end
p_hat = p_hat/n;	% normalize

%% Deblurred unimodal decomposition and JSD optimization
deblurred = unidec(p_hat,0);
assert(size(deblurred,1)==M,'bad M'); 
deblurred = deblurred(:,2:end-1);   % see plot shift below
deblurred = unimixmaxjsd(x,deblurred);

%% Reblur
reblurred = deblurred;
x0 = median(x);
blur = exp(-.5*((x-x0)/delta_h).^2)/(delta_h*sqrt(2*pi));
for m = 1:M
    reblurred(m,:) = conv(deblurred(m,:),blur,'same')*mean(diff(x));
end

%% Output (1)
% NB. Adding any additional fields here will break the plot macros below.
% TDE-relevant fields are added after the plot.
y.sweep = sweep;
y.unblurred = unblurred;
y.deblurred = deblurred;
y.reblurred = reblurred;

%% Plot
% Using macros and evals saves 60 lines of code. Also uses tight_subplot,
% available at https://www.mathworks.com/matlabcentral/fileexchange/27991
if plotflag
    M1 = max(M-1,1);
    padding = 1.1;
    maxStack = max([max(sum(sweep,1)),max(sum(unblurred,1)),...
        max(sum(deblurred,1)),max(sum(reblurred,1))])*padding;
    maxSuper = max([sweep(:);unblurred(:);deblurred(:);reblurred(:)])*padding;
    fields = fieldnames(y);
    %% Stacked components
    figure;
    ha = tight_subplot(numel(fields),1,.01,.1,.1);
    for j = 1:numel(fields)
        eval(['axes(ha(',num2str(j),'));']);
        eval(['a = area(x,',fields{j},''');']);
        eval(['for m = 1:M,','a(m).FaceColor = [(m-1)/M1,0,1-(m-1)/M1];',...
            'end,','set(gca,''XTick'',[],''YTick'',[]);',...
            'if plotflag > 0, legend(cellstr(num2str((1:M)''))''); end,',...
            'box on;','axis([min(x),max(x),0,maxStack])']);
        ylabel(fields{j},'Interpreter','latex');
    end
    %% Superimposed components
    figure;
    ha = tight_subplot(numel(fields),1,.01,.1,.1);
    for j = 1:numel(fields)
        eval(['axes(ha(',num2str(j),'));']);
        eval(['hold on;','for m = 1:M,','plot(x,',fields{j},...
            '(m,:),''Color'',[(m-1)/M1,0,1-(m-1)/M1],''LineWidth'',1);','end']);
        eval(['set(gca,''XTick'',[],''YTick'',[]);',...
        'if plotflag > 0, legend(cellstr(num2str((1:M)''))''); end,',...
        'box on;','axis([min(x),max(x),0,maxSuper])']);
        ylabel(fields{j},'Interpreter','latex');
    end
end

%% Output (2)
y.h = tde.h;
y.x = tde.x;
y.y = tde.y;
y.uc = tde.uc;
y.mfuc = tde.mfuc;
y.a = tde.a;
y.l = tde.l;
y.u = tde.u;