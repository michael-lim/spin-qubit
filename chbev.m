function ch = chbev(A)
    [~,n] = size(A);
    CU = lapack('zheev', 'V','U',n,A,n,zeros(n,1),zeros(2*n-1,1),2*n-1,zeros(3*n-2,1),0);
    [W,Z] = CU{[4,6]};
    D = zeros(n,n);
    for i=1:n
        D(i,i)=exp(Z(i));
    end
    ch = W*D/W;
end