function res = TestStability(A)
[~,J] = jordan(A); % computes the real Jordan canonical form

e = eig(J); % returns the eigenvalues of J;
n = length(A); % returns the length of the matrix A

if any(real(e) > 0)
    disp('The matrix A is unstable')
    res = 1;
elseif all(real(e) < 0)
    disp('The matrix A is asymptotically stable')
    res = -1;
else    
    % any(real(e) == 0);
    index = find(real(e)==0);
    m = length(index);    
    i = 1;
    while i <= m
        if (imag(e(index(i)))) == 0 % null eigenvalue
            if index(i)== n
                disp('The Matrix A is simply stable')
                res = 0;
                break
            elseif all(J(index(i):index(i)+1,index(i):index(i)+1) == [0 1 ; 0 0])
                disp('The Matrix A is unstable')
                res = 1;
                break
            end
            i = i+1;
        else % complex conjugate eigenvalues
            temp = J(index(i):index(i)+1,index(i):index(i)+1);
            if index(i) > n-3
                disp('The Matrix A is simply stable')
                res = 0;
                break
            elseif all(J(index(i):index(i)+3,index(i):index(i)+3) == [temp eye(2); zeros(2,2) temp])
                disp('The Matrix A is unstable')
                res = 1;
                break
            end
            i = i+2;
        end
    end
end
    
end