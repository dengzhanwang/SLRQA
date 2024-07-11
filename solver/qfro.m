function out=qfro(A)
out=0;
for ii=1:4
    B=A(:,:,ii);
    B=B(:);
    out=out+sum(B.*B);
end
end

