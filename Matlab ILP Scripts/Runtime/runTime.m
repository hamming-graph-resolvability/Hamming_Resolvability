tests = readTSV('res_set_data_3.tsv');
n = length(tests);
num_reps = 5;
for i = 1:n
    fprintf('Starting Example %i out of %i \n', i,n)
    times = zeros(1,num_reps);
    k = tests{i}{1};
    a = tests{i}{2};
    set = tests{i}{3};
    isResolving = tests{i}{4};
    for j = 1:num_reps
        tic
        A = makeMatrix(set,k,a);
        [resolve,x] = ILP_resolve(k,a,A);
        timeElapsed = toc;
        if isResolving ~= resolve
            fprintf('Incorrect answer: Correct = %i, Got = %i',isResolving,resolve)
            return
        end
        times(j) = timeElapsed;
    end
    tests{i}{5} = times;
end