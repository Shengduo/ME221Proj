%% Course project for HW1 in ME221
clc,clear;
close all;

% Define the geometry
Y = [-1/2, 1/2];
ns = [0.2, 0.05, 0.01];
mOverNs = [1., 2.];
NofPts = 100;
h = (Y(2) - Y(1)) / (NofPts - 1);

% Construct the grid
x1 = linspace(Y(1), Y(2), NofPts);
x2 = linspace(Y(1), Y(2), NofPts);
[X2, X1] = meshgrid(x2, x1);
alpha = ones(NofPts, NofPts);

%% Sub 1
% Loop through all cases
for i = 1:1:size(ns, 2)
    for j = 1:1:size(mOverNs, 2) 
        % Get m and n
        n = ns(i);
        m = mOverNs(j) * n;
        
        % Update alpha
        ellipseIdx = ((X1 ./ m) .^ 2 + (X2 ./ n) .^ 2 <= 1.);
        alpha(ellipseIdx) = 2.;
        
        % Loop through the xi's
        for d = 1:1:2
            % Compute div(alpha xi)
            divXi(d, :, :) = computeDivXi(alpha, h, d, NofPts);
        end
        % Solve the problem using fourier transform
        aStar = effectiveSolve(X1, Y1, alpha, NofPts);
    end
end

%% Function to get aStar
function aStar = effectiveSolve(X1, X2, alpha, NofPts)
    alpha0 = ones(NofPts, NofPts);

end

%% Function to compute divXi
function divXi = computeDivXi(alpha, h, d, NofPts)
    % a Xi = \alpha e_d
    % div(a Xi) = da / dx_d
    if d == 1 % X1 direction
        divXi = (alpha([2:1:NofPts, 1], :) - alpha([NofPts, 1:1:NofPts - 1], :)) ./ 2 / h;
    elseif d == 2 % X2 direction
        divXi = (alpha(:, [2:1:NofPts, 1]) - alpha(:, [NofPts, 1:1:NofPts - 1])) ./ 2 / h;
    end
end

%% Function to compute div(a \cdot grad(U))
function divAGradU = computeDivAGradU(alpha0, alpha, U, h, NofPts)
    % First Grad U
    GradU = zeros(2, NofPts, NofPts);
    
    % Compute the two directions
    divAGradU = zeros(NofPts, NofPts);
    aVals = alpha0 - alpha;
    
    % In X1 direction
    GradU(1, :, :) = (U([2:1:NofPts, 1], :) - U([NofPts, 1:1:NofPts - 1], :)) ./ 2 / h;

    % In X2 direction
    GradU(2, :, :) = (U(:, [2:1:NofPts, 1]) - U(:, [NofPts, 1:1:NofPts - 1])) ./ 2 / h;

    

end