function [C] = HRMC(X, G, alpha, beta, lambda, gamma, variance, tolX, maxiter)
    %%  HRMC High-rank matrix completion for gene prioritisation
    %   INPUT:
    %     - X: gene-disease association matrix.
    %     - G: graph containing disease-disease similarities.
    %     -- Parameters of the model: 
    %         -- alpha: weights the contribution of G.
    %         -- beta: l2 regularization penalty
    %         -- lambda: l1 regularization penalty.
    %         -- gamma: null-diagonal constraint penalty.
    %         -- variance: initial variance of the matrix.
    %    -- Stopping criteria:
    %        -- tolX: relative tolerance in the change in C.
    %        -- maxiter: maximun number of iterations.
    %  OUTPUT:
    %    - C: matrix of disease-disease associations.
    
    % initialisation
    [~, m] = size(X);
    C = rand(m, m) * sqrt(variance);
    C0 = C;
    I_W = eye(m); % identity matrix


    % iteration
    COV = X' * X;
    numerator = COV + alpha .* G;

    for iter = 1:maxiter
        % learn W
        epsilon = eps(C);
        denominator = (COV + alpha + beta) * C + lambda + gamma * I_W + epsilon;
        C = max(0, C .* (numerator ./ denominator));

        % get the max change in W
        dW = max(max(abs(C - C0) / (sqrt(eps) + max(max(abs(C0))))));
        C0 = C;

        % break point
        if dW <= tolX
            fprintf('\n iter %d, tolX %f\n', iter, full(dW));
            break;
        end


    end
end

