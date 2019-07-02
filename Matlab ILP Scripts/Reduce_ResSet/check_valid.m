function valid = check_valid(A,A_rank,X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
new_rank = rank(A);
if new_rank <= A_rank
    valid = false;
    return
end
if X
    b = A*X;
    if b == 0
        valid = false;
        return
    end
end
valid = true;
end

