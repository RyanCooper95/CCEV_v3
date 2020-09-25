function StabDetCheck(A,B,C,D)
n = length(A);
[Ak,~,~,~,n1,n2,n3,n4]=kalmcd(A,B,C,D); % Ultimate Kalman Decomposition
if n1 + n3 == 0
    disp('The system is fully observable');
else
    disp('The system is not fully observable');
    index = [1:1:n1 n1+n2+1:1:n-n4];
    res = TestStability(Ak(index,index));
    switch res
        case -1
            disp('The unobservable part is asymtptotically stable')
        case 0
            disp('The unobservable part is simply stable')
        case 1
            disp('The unobservable part is unstable')
    end
end

if n3 + n4 == 0
    disp('The system if fully reachable');
else
    disp('The system is not fully reachable');
    res = TestStability(Ak(n-n3+n4:n,n-n3+n4:n));
    switch res
        case -1
            disp('The unreacahble part is asymtptotically stable')
        case 0
            disp('The unreacahble part is simply stable')
        case 1
            disp('The unreacahble part is unstable')
    end
    
end


end