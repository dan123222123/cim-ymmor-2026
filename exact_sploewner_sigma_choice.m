%% Setup SISO system
n = 6; A = diag(-n:-1); B = (1:n)'; C = 1:n;
K = n; % data matrix size/half number of moments
%% ERA
sigma = Inf;
M = Numerics.sploewner.build_exact_moments(sigma,A,B,C,2*K);
[Db,Ds] = Numerics.sploewner.build_sploewner(sigma,M,M,M,K);
ew = realize_inorder(Db,Ds); ERA_err = norm(ew-diag(A))
%% realize system on a grid of shifts
N = 321;
% N = 1281;
x = linspace(-10,5,N); y = linspace(-7.5,7.5,N);
[X,Y] = meshgrid(x,y); G = X + 1i*Y;
SPLoewner_err = zeros(N,N);
parfor i=1:N
    for j=1:N
        sigma = G(i,j);
        M = Numerics.sploewner.build_exact_moments(sigma,A,B,C,2*K);
        [Db,Ds] = Numerics.sploewner.build_sploewner(sigma,M,M,M,K);
        ew = realize_inorder(Db,Ds); SPLoewner_err(i,j) = norm(ew-diag(A));
    end
end
%% Heatmap of Ratio between ERA and SPLoewner
ls_eravspl = log10(ERA_err./SPLoewner_err);
cl = max([ceil(max(ls_eravspl,[],"all")), abs(floor(min(ls_eravspl(~isinf(ls_eravspl)),[],"all")))]);
h = heatmap(x,y,ls_eravspl); clim([-cl cl]); colormap(redblue(5000));
% make better labels
CustomXLabels = string(x); CustomYLabels = string(y);
CustomXLabels(mod(x,1) ~= 0) = " "; CustomYLabels(mod(y,0.5) ~= 0) = " ";
h.XDisplayLabels = CustomXLabels; h.YDisplayLabels = CustomYLabels;
set(get(gca,'xlabel'),'rotation',90)
grid off;

function ew = realize_inorder(Db,Ds)
    ew = eig(Ds,Db); [~,ewidx] = sort(abs(ew),"descend"); ew = ew(ewidx);
end
