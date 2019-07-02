function A = makeMatrix(set,k,a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
set = strsplit(set,',');
nNodes = length(set);
A = zeros(nNodes, k*a);
for i = 1:nNodes
    str = set{i};
    encoding = zeros(a*k,1);
    for j = 1:k
        c = str2num(str(j))+1+a*(j-1);
        encoding(c) = 1;
        A(i,:) = encoding;
    end
end

end

