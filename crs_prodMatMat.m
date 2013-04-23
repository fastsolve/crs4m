function [C_rowptr,C_colind,C_val] = crs_prodMatMat(A_rowptr,A_colind,A_val,...
    B_rowptr,B_colind,B_val,N,M,L)
% Matrix-matrix multiplication C=A*B in compressed row format.
% Matrix A is NxM, Matrix B is MxL, and C is NxL.
% Algorithm is based on SMMP (http://www.netlib.org/aicm/smmp).

% %#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1]), coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1]),int32(0),int32(0),int32(0)}


%The whole process contains three steps
%1. Determine C_rowptr
%2. Fill in C_colind
%3. Fill in C_val

if nargin<9
    N=int32(length(A_rowptr))-1;
    M=int32(max(A_colind));
    L=int32(max(B_colind));
end

%initialization of row vector. It has the size of A_rowptr.
C_rowptr = nullcopy(zeros(length(A_rowptr),1,'int32'));
C_rowptr(1)=1;

index=zeros(max(M,L),1,'int32');

%Step one: Determine C_rowptr
for i=1:N
    istart=int32(-1);
    clength=int32(0);
    for jj=A_rowptr(i):A_rowptr(i+1)-1
        j=A_colind(jj);
        for k=B_rowptr(j):B_rowptr(j+1)-1
            if index(B_colind(k))==0
                index(B_colind(k))=istart;
                istart=B_colind(k);
                clength=clength+1;
            end
        end
    end
    C_rowptr(i+1)=C_rowptr(i)+clength;
    for j=C_rowptr(i):C_rowptr(i+1)-1
        k=istart;
        istart=index(istart);
        index(k)=0;
    end
    index(i)=0;
end

%step two: Fill in C_colind
C_colind = nullcopy(zeros(C_rowptr(N+1)-1,1,'int32'));

for i=1:N
    istart=int32(-1);
    clength=int32(0);
    for jj=A_rowptr(i):A_rowptr(i+1)-1
        j=A_colind(jj);
        for k=B_rowptr(j):B_rowptr(j+1)-1
            if index(B_colind(k))==0
                index(B_colind(k))=istart;
                istart=B_colind(k);
                clength=clength+1;
            end
        end
    end
    C_rowptr(i+1)=C_rowptr(i)+clength;
    for j=C_rowptr(i):C_rowptr(i+1)-1
        C_colind(j)=istart;
        istart=index(istart);
        index(C_colind(j))=0;
    end
    index(i)=0;
end

%The final step: Fill in values of C_val
C_val = nullcopy(zeros(C_rowptr(N+1)-1,1, class(A_val)));
temp  = zeros(length(index),1,class(A_val));

for i=1:N
    for jj=A_rowptr(i):A_rowptr(i+1)-1
        j=A_colind(jj);
        ajj=A_val(jj);
        for k=B_rowptr(j):B_rowptr(j+1)-1
            temp(B_colind(k))=temp(B_colind(k))+ajj*B_val(k);
        end
    end
    for j=C_rowptr(i):C_rowptr(i+1)-1
        C_val(j)=temp(C_colind(j));
        temp(C_colind(j))=0;
    end
end

function test  %#ok<DEFNU>
%!test
%! for k=1:100
%!     A = sprand(20,10,0.1);
%!     B = sprand(10,30,0.1);
%!     C = A*B;
%!     [A_rowptr,A_colind,A_val] = crs_createFromSparse(A);
%!     [B_rowptr,B_colind,B_val] = crs_createFromSparse(B);
%!     [C_rowptr, C_colind, C_val] = crs_prodMatMat(A_rowptr,A_colind,A_val, B_rowptr,B_colind,B_val,int32(20),int32(10),int32(30));
%!     assert( normest(C-crs_createSparse(C_rowptr, C_colind, C_val, 20, 30))<1.e-12);
%! end
