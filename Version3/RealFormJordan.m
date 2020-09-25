 function [V,J] = RealFormJordan(A)
[Vprime,~] = jordan(A); % computes the Jordan Canonical/Normal Form of the  matrix A.
% It also computes the similarity transformation, V, so that V\A*V = J.
% The columns of V are the generalized eigenvectors.
n = length(A); % returns the length of the matrix A
V = zeros(n,n); % is an n-by-n matrix of zeros

%% Look for indepedent eigenvectors
Vi = Vprime(:,1); % initialise the matrix of independent eigenvectors
i = 2; % start looking from the 2nd
while i <= n
    test = 0; % check if there is a new vector or not
    [~,m] = size(Vi); % number of columns in Vi
    for k = 1:m
        % check if the vectors are equal
        if all(real(Vprime(:,i)) == real(Vi(:,k))) && all(abs(imag(Vprime(:,i))) == abs(imag(Vi(:,k))))
            test = 1;
        end
    end
    if test == 0 % add the new vector to the matrix
        Vi = [Vi Vprime(:,i)];
    end
    i = i+1;
end
%% Reshafle the Eigenvectors Matrix
k = 1;
i = 1;
while i <= n
    if any(imag(Vi(:,k)))
        V(:,i) = real(Vi(:,k));
        V(:,i+1) = imag(Vi(:,k));
        i = i+2;
    else
        V(:,i) = Vi(:,k);
        i = i+1;
    end
    k=k+1;
end
%% Determine the Jordan Real Form (it works also in case of chains of generalised eigenvectors)
J = V\A*V;