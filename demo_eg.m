% This script is free to use and modify, and is openly distributed, 
% but the copyright is owned by Dr Stephen Suryasentana (stephen@autogeolab.com).
% This script may not be re-distributed as a part of a commercial product unless agreed upon with the copyright owner.
% This script is distributed WITHOUT ANY WARRANTY, or any implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% This script requires SeDuMi and YALMIP to be installed in Matlab.
%
% Any use of this script (or forks or modified versions of it) in a
% published work must cite the following publication:
% 'Suryasentana, S. K., Burd, H. J., Byrne, B. W., & Shonberg, A. (2020). A 
% systematic framework for formulating convex failure envelopes in multiple 
% loading dimensions. Géotechnique, 70(4), 343-353'.

clc;
close all;

% read in data
tbl = readtable('demo_eg.csv');

%% 1: Set up normalised data
F = [tbl.H tbl.M];

%% 2: Set up polynomial
[ndata, ndim] = size(F);
polydeg = 4; % degree of polynomial
x = sdpvar(ndim,1); % to contain the load variables e.g. H & M
b = sdpvar(ndata,1); % to contain data point evaluations 

% define polynomial in terms of H & M with unknown coefficients
[poly,a,monomials] = polynomial(x, polydeg, polydeg);

% display trial polynomial
sdisplay(poly)

%% 3: Set values to some of the known coefficients
cx = []; % to contain all constraints for optimisation
m1 = []; % to contain id of coefficients which has been set to 1
for i = 1:size(monomials, 1)
    degs = degree(monomials(i), x);

    % set coefficient to 1 for uniaxial capacity terms
    if any(degs == polydeg) 
        poly = replace(poly, a(i), 1);
        m1 = [m1;i];
    end    
end

%% 4: Apply constraints
% apply SOS constraint to Hessian of polynomial
H = hessian(poly, x);
na = size(a, 1);
monomials = sdisplay(monomials);
cx = [cx;sos(H);];

% apply data points to polynomial
for i = 1:ndata
    xi = F(i, :)'; % a single H and M dataset
    v2 = replace(poly, x, xi); % substitute H and M values into polynomial
    cx = [cx; abs(v2 - 1) <= b(i)]; % use a variable b as upper bound of polynomial evaluations of data point values
end

%% 5: Define objective function
obj = norm(b, 2);

%% 6: Solve optimisation problem by minimising objective function
[sol] = optimize(cx, obj,sdpsettings('solver','sedumi'))
fprintf('Minimised objective value = %.4f\n\n', value(obj));

%% 7: Display results
% display optimised values of coeffs (to 2 decimal places)
coeffs = round(value(a), 2);
coeffs(m1) = 1; % restore coeffs predefined as 1
coeffs

% display polynomial with optimised values of coeffs 
pf = replace(poly, a, coeffs);
sdisplay(pf)

%% 8: Plot results
% manually define yield surface formulation based on above results
yieldf = @(H, M) H.^4 -0.31*H.^3.*M + 0.8*H.^2.*M.^2 -1.39*H.*M.^3 + M.^4 - 1;

% plot yield surface and data points
scatter(F(:,1), F(:,2), 'o');
hold on
[x, y] = meshgrid( -1.4:0.01:1.4, -1.4:0.01:1.4);
z = yieldf(x, y);
contour(x,y,z, [0, 0]);