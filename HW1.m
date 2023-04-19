%% Course project for HW1 in ME221
clc, clear;
close all;

% Define the geometry
Y = [-1/2, 1/2];

% ------------------- For sub problem 1 ------------------- 
ns = [0.2, 0.05, 0.01];
mOverNs = [1., 2.];
% --------------------------------------------------------- 

% ------------------- For sub problem 2 ------------------- 
ns = [50];
mOverNs = [0.005];
% --------------------------------------------------------- 

NofPts = 1024;
h = (Y(2) - Y(1)) / (NofPts);

% Construct the grid
x1 = linspace(Y(1), Y(2), NofPts + 1);
x2 = linspace(Y(1), Y(2), NofPts + 1);
x1 = x1(1 : end - 1);
x2 = x2(1 : end - 1);
[X2, X1] = meshgrid(x2, x1);
alpha = ones(NofPts, NofPts);

% Parameters to start the solve
U0 = ones(NofPts, NofPts);
alpha0 = 1.;
rtol = 1.e-15;
max_iters = 10000;
nOfCases = length(ns) * length(mOverNs);
aStars = zeros(nOfCases, 2, 2);
nns = zeros(1, nOfCases);
mms = zeros(1, nOfCases);

%% Sub 1
% Loop through all cases
caseNo = 1;
for i = 1:1:size(ns, 2)
    for j = 1:1:size(mOverNs, 2) 
        % Get m and n
        n = ns(i);
        m = mOverNs(j) * n;
        nns(caseNo) = n;
        mms(caseNo) = m;

        % Update alpha
        alpha = ones(NofPts, NofPts);
        ellipseIdx = ((X1 ./ m) .^ 2 + (X2 ./ n) .^ 2 <= 1.);
        alpha(ellipseIdx) = 2.;
        spy(alpha - 1);
        divXis = zeros(2, NofPts, NofPts);
        r_errs = zeros(1, 2);
        Us = zeros(2, NofPts, NofPts);
        ress = zeros(1, 2);

        % Loop through the xi's
        for d = 1:1:2
            % Compute div(alpha xi)
            divXi = computeDivXi(alpha, h, d, NofPts);
            divXis(d, :, :) = divXi;
            [r_err, U] = solveU(alpha0, alpha, NofPts, h, divXi, ...
                                U0, rtol, max_iters);
            r_errs(d) = r_err;
            Us(d, :, :) = U;
            res = computeRes(alpha, U, h, NofPts, divXi);
            disp(strcat("Residual = ", num2str(res)));
        end
        
        % Compute aStar
        aStar =  zeros(2, 2);
        for ii = 1:1:2
            for jj = 1:1:2
                U = squeeze(Us(jj, :, :));
                if ii == 1
                    Ujk = (U([2:NofPts, 1], :) - U([NofPts, 1:NofPts - 1], :)) ./ (2*h);
                    AikUjk = alpha .* Ujk;
                elseif ii == 2
                    Ujk = (U(:, [2:NofPts, 1]) - U(:, [NofPts, 1:NofPts - 1])) ./ (2*h);
                    AikUjk = alpha .* Ujk;
                end
                aStar(ii, jj) = sum(((ii == jj) .* alpha + AikUjk) .* h^2, 'all');
            end
        end
        disp("-------------------------------------------------");
        disp(strcat("Case No: ", num2str(caseNo)))
        disp(strcat("m, n = ", num2str([m, n])));
        disp("aStar:");
        disp(aStar)
        aStars(caseNo, :, :) = aStar;
        caseNo = caseNo + 1;

        % Solve the problem using to get Us
    
    end
end


%% Solve for U iteratively
function [r_err, U] = solveU(alpha0, alpha, NofPts, h, divXi, U0, rtol, max_iters)
    % Initialize U
    U = U0;
    r_err = 1.;
    iter = 0;
    
    % The wave number matrix
    L = h * (NofPts);
    % K = zeros(NofPts, NofPts);
    
%     for i = 1:1:NofPts
%         for j = 1:1:NofPts
%             K(i, j) = -((i - 1) ^ 2 + (j - 1) ^ 2) / (L ^ 2) * 4 * pi^2;
%         end
%     end
    % k = (2*pi/L) * [0:1:NofPts - 1];
    k = (2*pi/L)*[0:(NofPts/2-1) (-NofPts/2):(-1)]; % Vector of wavenumbers
    [KX, KY]  = meshgrid(k,k); % Matrix of (x,y) wavenumbers corresponding
                                % to Fourier mode (m,n)
    delsq = -(KX.^2 + KY.^2); % Laplacian matrix acting on the wavenumbers
    delsq(1,1) = 1;           % Kluge to avoid division by zero of (0,0) waveno.
                                % of fhat (this waveno. should be zero anyway!)
                            
    % K = 1. ./ delsq;
    % K(1, 1) = 0.;
    while (iter < max_iters) && (r_err >= rtol)
        % Calculate f
        f = computeF(alpha0, alpha, U, h, NofPts, divXi);
        fHat = fft2(f ./ (-alpha0));
        
        % Solve for U;
        U_prev = U;
        U = ifft2(fHat ./ delsq);
        U = real(U);
        U = U - U(1, 1);
        
        iter = iter + 1;
        r_err = norm(U - U_prev, 'fro') / norm(U_prev, 'fro');
        
%         % Check residual now
%         Uxx = (U([2:NofPts, 1], :) - 2 .* U + U([NofPts, 1:NofPts - 1], :)) / (h ^ 2);
%         Uyy = (U(:, [2:NofPts, 1]) - 2 .* U + U(:, [NofPts, 1:NofPts - 1])) / (h ^ 2);
%         res = (Uxx + Uyy) + f ./ alpha0;
%         disp(strcat("Res in solver: ", num2str(sum(res.^2 .* h^2, 'all'))));
    end
    
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
    % Compute the two directions
    aVals = alpha0 - alpha;
    
    % In X1 direction
    AGradU(1, :, :) = aVals .* (U([2:1:NofPts, 1], :) - U([NofPts, 1:1:NofPts - 1], :)) ./ 2 / h;

    % In X2 direction
    AGradU(2, :, :) = aVals .* (U(:, [2:1:NofPts, 1]) - U(:, [NofPts, 1:1:NofPts - 1])) ./ 2 / h;
    
    divAGradU = (AGradU(1, [2:1:NofPts, 1], :) - AGradU(1, [NofPts, 1:1:NofPts - 1], :)) ./ (2 * h) ...
                + (AGradU(2, :, [2:1:NofPts, 1]) - AGradU(2, :, [NofPts, 1:1:NofPts - 1])) ./ (2 * h);
    
    divAGradU = squeeze(divAGradU);
end


%% Function to compute f
function f = computeF(alpha0, alpha, U, h, NofPts, divXi)
    f = -computeDivAGradU(alpha0, alpha, U, h, NofPts) + divXi;
end

%% Function to compute residual
function res = computeRes(alpha, U, h, NofPts, divXi)
    % First calculate divAlphaGradU
    aVals = alpha;
    AGradU = zeros(2, NofPts, NofPts);
    % In X1 direction
    AGradU(1, :, :) = aVals .* (U([2:1:NofPts, 1], :) - U([NofPts, 1:1:NofPts - 1], :)) ./ 2 / h;

    % In X2 direction
    AGradU(2, :, :) = aVals .* (U(:, [2:1:NofPts, 1]) - U(:, [NofPts, 1:1:NofPts - 1])) ./ 2 / h;
    
    divAGradU = (AGradU(1, [2:1:NofPts, 1], :) - AGradU(1, [NofPts, 1:1:NofPts - 1], :)) ./ (2 * h) ...
                + (AGradU(2, :, [2:1:NofPts, 1]) - AGradU(2, :, [NofPts, 1:1:NofPts - 1])) ./ (2 * h);
    
    divAGradU = squeeze(divAGradU);
    
    % Return residual
    res = sum((divAGradU + divXi) .^ 2 .* h^2, 'all');
end