/*
 * Some kernel macros for sparse matrices
 */

#ifndef _MSP_KERNEL_
#define _MSP_KERNEL_

/* Computes t = x[0:n-1]'*y[i[0:n-1]-1] */
__inline__ double MSP_ddot(const double * restrict x, 
                           const int32_T * restrict i, 
                           const double * restrict y, 
                           int32_T n) {
    int32_T j;
    double t = 0.0;
    
    y=y-1;
    for (j=0; j<n; ++j) t += (*x++) * y[*i++];
    return t;
}

/* Computes t = x[0:n-1]'*y[i[0:n-1]-1] */
__inline__ float MSP_sdot(const float * restrict x,
                          const int32_T * restrict i,
                          const float * restrict y,
                          int32_T n) {
    int32_T j;
    float t = 0.0;
    
    y=y-1;
    for (j=0; j<n; ++j) t += (*x++) * y[*i++];
    return t;
}

/* Computes y[i[0:n-1]-1] += a*x[0:n-1] */
__inline__ void MSP_saxpy(float a, 
                          const float * restrict x, 
                          float * restrict y, 
                          const int32_T * restrict i, 
                          int32_T n) {
    int32_T j=1;
    y = y-1;
    for (j=0; j<n; ++j) y[*i++] += a * (*x++);
}

/* Computes y[i[0:n-1]-1] += a*x[0:n-1] */
__inline__ void MSP_daxpy(double a, 
                          const double *restrict x, 
                          double *restrict y, 
                          const int32_T *restrict i, 
                          int32_T n) {
    int32_T j=1;
    y = y-1;
    for (j=0; j<n; ++j) y[*i++] += a * (*x++);
}

#endif
