function [output, isCRS] = CuBCRS(varargin) %#codegen
%Createa an opaque sparse matrix object in block-CRS format on CUDA.
%
%  CuBCRS() returns a definition of a structure, suitable in the
%  argument specification for codegen. (Note: Any compiled code
%  with CuBCRS as an input variable in its top-level would support
%  only double precision.)
%
%  CuBCRS(rowptr, colind, vals, [msp_type], [mb], [nb], [nnzb], [blkdim]) 
%  wraps the given opaque CUDA pointers into sparse mb-by-nb block-matrix 
%  in CuBCRS format. Each square block is of size blkdim-by-blkdim, stored 
%  in column major. If type is not given, it inherits from vals. If mb is 
%  not given, it is set to length(rowptr)-1. If nb is not given, it is 
%  assumed to be the same as mb. If nnzb is not given, it is set to 
%  length(colind). If blkdim is not given, then it is assumed to be 1, 
%  and the data structure is then essentially the same as CuCRS.
%
%  CuBCRS(obj, 'rowptr') or CuBCRS(obj, 'colind') converts the rowptr or
%  colind field into a CUDA pointer, respectively.
%
%  CuBCRS(obj, 'vals', type) converts the vals field into a CUDA pointer
%  of a particular pointer type specified by type, such as 'double *'.
%
% See also CuCRS, CuVec, cuBCRSCreate

coder.inline('always');

narginchk(0, 8);

% Use the second output argument to control customization for CRS
isCRS = nargout==2;
if isCRS
    dimsb_name = 'dims';
    nnzb_name = 'nnz';
    cstructname = 'MSP_CuCRS';
else
    dimsb_name = 'dimsb';
    nnzb_name = 'nnzb';
    cstructname = 'MSP_CuBCRS';
end

if nargin==0 
    if isCRS
        output = coder.typeof(struct(...
            'rowptr', uint64(0), 'colind', uint64(0), ...
            'vals', uint64(0), 'type', int32(0), ...
            dimsb_name, m2c_dims(2), nnzb_name, int32(0)));
    else
        output = coder.typeof(struct(...
            'rowptr', uint64(0), 'colind', uint64(0), ...
            'vals', uint64(0), 'type', int32(0), ...
            dimsb_name, m2c_dims(2), nnzb_name, int32(0), ...
            'blkdim', int32(0)));
    end
    coder.cstructname(output, cstructname);
elseif nargin==2 && ~isequal(varargin{2}, 'vals')
    output = m2c_castdata('int *', varargin{1}.(varargin{2}));
elseif nargin==3 && isequal(varargin{2}, 'vals')
    output = m2c_castdata(varargin{3}, varargin{1}.(varargin{2}));
elseif nargin>=3
    output = struct(...
        'rowptr', uint64(0), 'colind', uint64(0), ...
        'vals', uint64(0), 'type', int32(0), ...
        dimsb_name, [int32(0), int32(0)], nnzb_name, int32(0));
    if ~isCRS; output.blkdim = int32(1); end
    coder.cstructname(output, cstructname);
    
    output.rowptr = varargin{1}.data;
    output.colind = varargin{2}.data;
    output.vals = varargin{3}.data;
    
    if nargin>=4
        output.type = int32(varargin{4});
    else
        output.type = varargin{3}.data;
    end
    
    if nargin>=5
        output.(dimsb_name)(1) = int32(varargin{5});
    else
        output.(dimsb_name)(1) = varargin{1}.len-1;
    end
    
    if nargin>=6
        output.(dimsb_name)(2) = int32(varargin{6});
    else
        output.(dimsb_name)(2) = output.(dimsb_name)(1);
    end
    
    if nargin>=7
        output.(nnzb_name) = int32(varargin{7});
    else
        output.(nnzb_name) = int32(varargin{2}.len);
    end
    
    if isfield(output, 'blkdim') && nargin>=8
        output.blkdim = int32(varargin{8});
    end
else
    % Undefined.
end
