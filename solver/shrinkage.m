function output=shrinkage(A,epsilon)
    sign_A=sign(A);
    output=sign_A.*max(abs(A)-epsilon,0);
end