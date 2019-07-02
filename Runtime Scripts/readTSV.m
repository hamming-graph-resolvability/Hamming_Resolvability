function tests = readTSV(file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
data = tdfread(file);
n = length(data.k);
tests = {};
counter = 0;
for i = 1:n
    k = data.k(i);
    a = data.a(i);
    isResolving = strtrim(data.resolving(i,:));
    if strcmp(isResolving,'True')
        isResolving = true;
    else
        isResolving = false;
    end
    set = strtrim(data.set(i,:));
    if a*k <= 25
        counter = counter+1;
        tests{counter} = {k,a,set,isResolving,0};
    end
end
end

