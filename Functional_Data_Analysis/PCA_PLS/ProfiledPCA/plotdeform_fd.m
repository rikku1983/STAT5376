function Hline = plotdeform_fd(Wfdobj, Lfdobj, matplt, href, nfine)
%  PLOTDEFORM   Plot a deformation function 
%            d(t) = h(t) - t
%  defined by functional data object Wfd, the log-derivative of
%  one or more  strictly monotonic time-warping functions 
%        h(t) = int_0^t exp[W(v)] dv/int_0^T exp[W(v)] dv
%  over interval [0,T]
%
%  Arguments:
%  FDOBJ   ... A functional data object to be plotted.
%  LFDOBJ  ... A linear differential operator object or a positive
%              integer specifying a derivative.  This operator is
%              applied to the functions before plotting.
%  MATPLOT ... If MATPLT is nonzero, all curves are plotted in a 
%              single plot.
%              Otherwise, each curve is plotted separately, and the
%              next curve is plotted when the mouse is clicked.
%  HREF    ... If HREF is nonzero, a horizontal dotted line through 0 
%              is plotted.
%  NX      ... The number of plotting points to be used.

%  Last modified 25 February 2012

%  set default arguments

if nargin < 4 || isempty(href),   href = 1;             end
if nargin < 3 || isempty(matplt), matplt = 1;           end
if nargin < 2 || isempty(Lfdobj), Lfdobj = int2Lfd(0);  end

%  check arguments

if ~isa_fd(Wfdobj)
    error ('Argument Wfdobj not a functional data object.');
end

Lfdobj = int2Lfd(Lfdobj);
if ~isa_Lfd(Lfdobj)
    error ('Argument Lfdobj not a linear differential operator object.');
end

% set up dimensions of problem

coef   = getcoef(Wfdobj);
coefd  = size(coef);
ndim   = length(coefd);
ncurve = coefd(2);

if ndim > 2
    error('Plotting is not available for multivariate functions.');
end

%  extract basis information

Wbasis  = getbasis(Wfdobj);
nWbasis = getnbasis(Wbasis);
rangex  = getbasisrange(Wbasis);
width   = diff(rangex);

%  set up fine mesh of evaluation points and evaluate curves

if nargin < 5
    nx = max([101, 10*nWbasis+1]);
end

x        = linspace(rangex(1),rangex(2),nx)';
fx       = monfn(x, Wfdobj);
fmax     = fx(nx,:);
hx       = x(1) + width.*fx./(ones(nx,1)*fmax);
hx(1,:)  = x(1);
hx(nx,:) = x(nx);
dx = hx - x*ones(1,ncurve);

%  calculate range of values to plot

frng(1) = min(min(dx));
frng(2) = max(max(dx));

%  fix range if limits are equal

if frng(1) == frng(2)
    if abs(frng(1)) < 1e-1
        frng(1) = frng(1) - 0.05;
        frng(2) = frng(2) + 0.05;
    else
        frng(1) = frng(1) - 0.05*abs(frng(1));
        frng(2) = frng(2) + 0.05*abs(frng(1));
    end
end

%  extract argument, case and variable names

fdnames  = getnames(Wfdobj);

if matplt
    for icurve=1:ncurve
        if icurve > 1
            hold on
        end
        dxi = dx(:,icurve);
        %  plot all curves
        if href && (frng(1) <= 0 && frng(2) >= 0)
            if nargout > 0
                Hline = plot (x, dxi, '-', x, zeros(nx), ':');
            else
                plot (x, dxi, '-', x, zeros(nx), ':')
            end
        else
            if nargout > 0
                Hline = plot (x, dx(:,icurve), '-');
            else
                plot (x, dx(:,icurve), '-')
            end
        end
        xlabel(['\fontsize{12} ',fdnames{1}]);
        if iscell(fdnames{3})
            ylabel(['\fontsize{12} ',fdnames{3}{1}])
        else
            ylabel(['\fontsize{12} ',fdnames{3}])
        end
        if frng(2) > frng(1)
            axis([x(1), x(nx), frng(1), frng(2)]);
        end
    end
    if ncurve > 1
        hold off
    end
else
    %  page through curves one by one
    for icurve = 1:ncurve
        if href && (frng(1) <= 0 && frng(2) >= 0)
            if nargout > 0
                Hline = plot(x, dx(:,icurve), 'b-', ...
                    [min(x),max(x)], [0,0], 'r:');
            else
                plot(x, dx(:,icurve), 'b-', ...
                    [min(x),max(x)], [0,0], 'r:')
            end
        else
            if nargout > 0
                Hline = plot(x, dx(:,icurve), 'b-');
            else
                plot(x, dx(:,icurve), 'b-')
            end
        end
        xlabel(['\fontsize{12} ',fdnames{1}])
        if iscell(fdnames{3})
            ylabel(['\fontsize{12} ',fdnames{3}{1}])
        else
            ylabel(['\fontsize{12} ',fdnames{3}])
        end
        if iscell(fdnames{2})
            title( ['\fontsize{12}', fdnames{2}{2}(icurve,:)]);
        else
            title( ['\fontsize{12} Curve ', num2str(icurve)]);
        end
        pause;
    end
end
