%Number of unknowns, M points in (x) and N points in (y)
function [u] = test(fPassed,MPassed,NPassed,cPassed)
M = MPassed;
N = NPassed;
size = (N-1)*(M-2);
c = cPassed;
dt = 1; %TODO: Calculate dt from the given paramater
dx = 3; % TODO: Calculate dx from the given parameter
A = 2*c*dt;
B = -4*c*dt;
C = dx^2;
D = -dx^2;  

%Matrix = zeros(size);
Matrix = sparse(size,size);

%Matrix of Coefficients 
for i = 1:size
    Matrix(i,i) = B;
        if i+(N-1) <= size
            Matrix(i,i+(N-1)) = A;
        end
        if i-(N-1) >= 1
            Matrix(i,i-(N-1)) = A;
        end
        if i >= 2
            %Coefficients of (i,j-1)
            if mod(i-1,(N-1)) == 0
                Matrix(i,i-1) = 0;
            else
                Matrix(i,i-1) = C;
            end
            
        end
        if i<= size -1
            %Coefficients of (i,j+1)
            if mod(i,(N-1)) == 0
                Matrix(i,i+1) = 0;
            else
                Matrix(i,i+1) = D;
            end
            
        end
end

% RHS and Calculating the solution
U = zeros (1,M-2);
for idx = 1 : M-2
    f = fPassed(dx *idx);
    U(1,idx) = f;
end
%fprintf('%g ', U);
len = (M-2) * (N-1);
b = zeros(len,1);
for idxM = 1:M-2
 
    if (idxM==1) 
    b(idxM,1) = U(1,idxM);
    else
        b(((idxM-1)*(N-1)+1),1) = U(1,idxM);
    end
end

b = -b;
b2 = sparse(b);
S = sparse(Matrix);
sol = S\b2;
%sol = Matrix\b;
%clear b Matrix c dt f dx A B C D c size len idxM idx i 
clear b Matrix c dt f dx A B C D c size len b2 S

% Concatonation Region to return the Matrix of the general solution
matrix = zeros(N,M);
submatrix = zeros(N-1,M-2);
for i=2:M-1
    matrix(N,i)= U(1,i-1);
end
sol = flip(sol);
cursor =1;
for i=1:M-2
   subcolumn = sol(cursor:(cursor-1)+N-1);
   cursor = cursor+N-1;
   submatrix(:,M-i)= subcolumn;
end
submatrix(:,1)=[];
matrix(1:N-1, 2:M-1)= submatrix;
u = matrix;
iteration =0;
while (iteration <2000)
    iteration = iteration +1;
    if( mod(iteration,10)==0;
        imgesc(u)
        colorbar;
        drawnow;
    end
    Lu=del2(u);
    u(2:N-1,2-M-1) = u (2:N-1,2:M-1) Lu(2:N-1,2:M-1);
clear M N submatrix cursor sol U subcolumn i matrix
end