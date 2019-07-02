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

%% Binary Search
top = 82;
bot = 1;
n = 800;
while top-bot > 1
    diff = top-bot;
    curr_size = bot+ceil(diff/2);
    fprintf('top = %i, bot = %i, curr_size = %i \n',top,bot,curr_size)
    tic
    found = false;
    for j = 1:n
        new_set = randperm(82,curr_size);
        A = A_82(new_set,:);
        [is_resolving,X] = ILP_resolve(k,a,A);
        if is_resolving
            set = cell(curr_size,1);
            for i = 1:curr_size
                set{i} = str{new_set(i)};
            end
            found = true;
            break
        end
    end
    toc
    if found
        top = curr_size;
    else
        bot = curr_size;
    end
end
fprintf('Last check for bot = %i\n',bot)
tic
for j = 1:n
        new_set = randperm(82,curr_size);
        A = A_82(new_set,:);
        [is_resolving,X] = ILP_resolve(k,a,A);
        if is_resolving
            fprintf('Smallest size is %i\n',bot)
            set = cell(curr_size,1);
            for i = 1:curr_size
                set{i} = str{new_set(i)};
            end
            break
        end
end  
toc