function [B_rowptr, B_colind, B_val] = crs_prodPtAP( A_rowptr, A_colind, A_val, P_rowptr, P_colind, P_val)
% Compute B=P'*A*P, where A and P are stored in compressed row storage.

% %#codegen -args {coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1]), coder.typeof(int32(0), [inf,1]), coder.typeof(int32(0), [inf,1]),
% %#codegen coder.typeof(0, [inf,1])}

% Computes B=P'*A*B;
[C_rowptr, C_colind, C_val] = crs_prodMatMat( A_rowptr, A_colind, A_val, P_rowptr, P_colind, P_val);
[Pt_rowptr, Pt_colind, Pt_val] = crs_transpose(P_rowptr, P_colind, P_val);
[B_rowptr, B_colind, B_val] = crs_prodMatMat( Pt_rowptr, Pt_colind, Pt_val, C_rowptr, C_colind, C_val);
