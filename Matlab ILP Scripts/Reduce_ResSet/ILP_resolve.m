function [resolving,X] = ILP_resolve(k,a,A)
%ILP_resolve takes k,a,A and checks if the set that formed A is resolving
%or not. If it finds a valid certificate that the set is NOT resolving it
%will return the valid solution.

    %fprintf('\n=== RUNNING ILP TO CHECK IF SET IS RESOLVING ===\n');
    d   = a*k;
    mat = @(x) reshape( x, a, k ); % turns (a*dk x 1 into a x k )

    % Pick a strategy
    if d <= 16
        c   = 2.^(1:d);  % for d > 16, we have numerical issues with this...
    else
        % Probability 1
        c = randn(1,d);
    end

    % if x has any non-zero entries, in range {-1,0,1},
    %   then dot( c, x ) must be non-zero because there
    %   are no cancellations.
    % Thus, if we have a non-zero feasible point x,
    %   then either dot(c,x) is > 0 or it is < 0
    % So, we have two cases
    %   (We have to do it this way since maximize sum(abs(x))
    %                           or even  maximize dot(c,x)
    %   isn't convex).

    %tic
    cvx_begin quiet
        cvx_solver gurobi % handles integer constraints
        variable x(d) integer
        %minimize c*x
        minimize 0
        subject to
            abs(x) <= 1 % so {-1,0,1} range
            sum( mat(x) ) == 0 % Want blocks to sum to 0
            sum( mat(abs(x) ) ) <= 2
            A*x == 0
            c*x <= -1e-3 % quickly find certificate
    cvx_end
    %toc
    if strcmpi( cvx_status, 'infeasible' )
        resolving = true;
    else
        resolving = false;
        if norm(x) < 1e-8
            fprintf('numerical errors')
        end
    end

    % If we claimed to be resolving, do a final check
%     if resolving
%         %disp('Feasiblity problem suggests this is a resolving set; running slower ILP to double-check');
%         if d > 16
%             % Draw a new random variable
%             c = randn(1,d);
%         end
%         %tic
%         cvx_begin quiet
%             cvx_solver gurobi % handles integer constraints
%             variable x(d) integer
%             minimize c*x
%             subject to
%                 abs(x) <= 1 % so {-1,0,1} range
%                 sum( mat(x) ) == 0 % Want blocks to sum to 0
%                 sum( mat(abs(x) ) ) <= 2
%                 A*x == 0
%         cvx_end
%         %toc
%         if ~all( x== 0 )
%             %disp('Slower ILP found a non-zero solution, so this is actually not resolving');
%             resolving = false;
%         end
%     end


    if ~resolving
        %fprintf(2,'Found a certificate that this is NOT a resolving set\n');
        X   = x;
    else
        %fprintf(2,'Only solution is the zero solution, so this IS a resolving set\n');
        X   = x;
    end

end

