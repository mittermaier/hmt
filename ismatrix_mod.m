function aus=ismatrix_mod(A)

%This function checks if the variable 'A' is a matrix. If so it returns 1
%If not it returns 0
%vectors are not deemed as matrices


if (ndims(A)==2)
    if (size(A)>ones(1,2))
        aus=1;
    else
        aus=0;
    end
else
    aus=0;
end
