close all;
clc;
clear all;


%% ReadMe
% To be completed later
% Improve the code, adding Polynomial Code and PolyDot

%% Input Parameters

% computational inputs, two N × N matrices A and B
N = 10; % The dimension (row/column) of each of the square matrices being multiplied.
A = randn(N,N);
B = randn(N,N);

m = 2; % The storage parameter 
% A fixed 1/m fraction of each input can be stored at each node.

matdot_k = 2 * m -1; % Matdot recovery threshold

P = 4; % The total number of worker nodes used.

%% Matdot Coded Computation
 
% Matrix A is split vertically and B is split horizontally
% A = [A_0 A_1 ... A_{m-1}]
% B = [B0 ; B1 ; .... ; B_{m-1}]

% Master Node Split
A_split = reshape ( A , N , N / m , m );
B_split = permute ( reshape (permute(B,[2 1]), N, N/m, m), [2 1 3]);

fusion_node_collect = zeros (N, N, P);
workers_flag = zeros (1 , P); % successful workers

Evaluation_points = randn(1,P);

for w = 1:1:P

    % Master node to w'th worker
    x = Evaluation_points(1,w);
    P_A = zeros ( N , N / m );
    P_B = zeros ( N / m , N );
    for i = 0:1:(m-1)
        A_i = zeros(N, N/m); A_i(:,:) = A_split (:,:,i+1);
        B_i = zeros(N/m, N); B_i(:,:) = B_split (:,:,i+1);
        P_A = P_A + ( x^i * A_i );
        P_B = P_B + ( x^(m-1-i) * B_i );
    end
    
    % No straggler for now
    % w'th worker computation
    fusion_node_collect (:,:,w) = P_A * P_B;
    workers_flag (1,w) = 1;
end

% Fusion node computation

% choosing k successful worker
k = matdot_k;
success_num = sum (workers_flag);
success_ind = find (workers_flag);
success_ind = success_ind ( randperm(success_num) );
output = zeros ( N , N );
if ( success_num >= k )
    for u = 1:1:N
        for v = 1:1:N
            Coeff = repmat(Evaluation_points(1,success_ind(1:k))',1,k).^ repmat(0:k-1,k,1);
            Eval = reshape ( fusion_node_collect (u,v,success_ind(1:k)) , k , 1);
            polynomial_coeff = Coeff\Eval;
            output ( u , v ) = polynomial_coeff(m);
        end
    end
    
else
    disp('Too many stragglers!')
end

sum(sum(output-(A*B)))
