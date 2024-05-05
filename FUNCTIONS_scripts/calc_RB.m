%compute the Algebraic Robustness

function algebraicConnectivity = calc_RB(A)
    % Calculate the Laplacian matrix
    degreeMatrix = diag(sum(A));
    laplacianMatrix = degreeMatrix - A;
    
    % Calculate the eigenvalues of the Laplacian matrix
    eigenvalues = eig(laplacianMatrix);
    
    % Sort the eigenvalues in ascending order
    sortedEigenvalues = sort(eigenvalues);
    
    % Find the second-smallest eigenvalue (excluding zero)
    algebraicConnectivity = sortedEigenvalues(2);
end
    