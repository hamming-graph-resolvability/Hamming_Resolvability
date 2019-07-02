%% Load the known length 82 resolving set for H_{8,20}
fid     = fopen('octamers.txt','r');
str     = textscan(fid,'%s');
str     = str{1};
fclose(fid);

nNodes = length(str);
k   = 8;
a   = 20;
A   = zeros( nNodes, k*a );
alphabet = unique( strcat( str{:} ) );
if length( alphabet ) ~= a, disp('Error, too many alphabet entries'); end
% Make a lookup for the one-hot encoding
alphabetNumeric = double( alphabet );
encoding        = zeros( max(alphabetNumeric), 1 );
encoding( alphabetNumeric ) = 1:a;

for i= 1:nNodes
    s   = str{i};
    z   = zeros(a,k);
    for j = 1:k
        c = s(j)+0;
        z( encoding(c), j ) = 1;
    end
    A(i,:) = z(:)';
end
A_82 = A;
tic
is_resolving = false;
set = cell(78,1);
counter = 0;
while ~is_resolving
    counter = counter+1;
    fprintf('Set %i\n',counter)
    inds = randperm(82,78);
    A = A_82(inds,:);
    [is_resolving,X] = ILP_resolve(k,a,A);
end
toc
fprintf('Finished')
for i = 1:78
    set{i} = str{inds(i)};
end
