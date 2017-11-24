function aus=iscolumn_mod(v)

%This function checks if the variable 'v' is a column vector. If so it returns 1
%If not it returns 0


if(isvector(v)==1)
    k=size(v);
    if (k(2)==1)
        aus=1;
    else
        aus=0;
    end
else
    aus=0;
end
