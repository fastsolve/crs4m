function flag = m2c_mkl
%Flag indicating whether m2c_mkl is on.
%It is always false within MATLAB. For m2c, it can be turned on with the
%-mkl option. It can also be turned on or off by the compiler flag
%-DM2C_MKL=1  DM2C_MKL=0, respectively.

coder.inline('always');
coder.cinclude('mspack.h');

if isempty(coder.target)
    flag = false;
else
    flag = int32(1); %#ok<NASGU>
    flag = coder.ceval(' ', coder.opaque('int', 'M2C_MKL'));
end
end
