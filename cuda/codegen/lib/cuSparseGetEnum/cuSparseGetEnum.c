#include "cuSparseGetEnum.h"
#include "mspack.h"
#include "m2c.h"

int cuSparseGetEnum(const emxArray_char_T *str)
{
  int val;
  boolean_T b_bool;
  int kstr;
  int exitg11;
  static const char cv0[26] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'P',
    'O', 'I', 'N', 'T', 'E', 'R', '_', 'M', 'O', 'D', 'E', '_', 'H', 'O', 'S',
    'T' };

  int exitg10;
  static const char cv1[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'P',
    'O', 'I', 'N', 'T', 'E', 'R', '_', 'M', 'O', 'D', 'E', '_', 'D', 'E', 'V',
    'I', 'C', 'E' };

  int exitg9;
  static const char cv2[23] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S' };

  int exitg8;
  static const char cv3[31] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T', 'I',
    'A', 'L', 'I', 'Z', 'E', 'D' };

  int exitg7;
  static const char cv4[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A', 'I',
    'L', 'E', 'D' };

  int exitg6;
  static const char cv5[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_', 'V',
    'A', 'L', 'U', 'E' };

  int exitg5;
  static const char cv6[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S', 'M',
    'A', 'T', 'C', 'H' };

  int exitg4;
  static const char cv7[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_', 'E',
    'R', 'R', 'O', 'R' };

  int exitg3;
  static const char cv8[32] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O', 'N',
    '_', 'F', 'A', 'I', 'L', 'E', 'D' };

  int exitg2;
  static const char cv9[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_', 'S',
    'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L', '_',
    'E', 'R', 'R', 'O', 'R' };

  int exitg1;
  static const char cv10[41] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'T', 'R', 'I', 'X', '_', 'T',
    'Y', 'P', 'E', '_', 'N', 'O', 'T', '_', 'S', 'U', 'P', 'P', 'O', 'R', 'T',
    'E', 'D' };

  b_bool = false;
  if (str->size[1] != 26) {
  } else {
    kstr = 0;
    do {
      exitg11 = 0;
      if (kstr + 1 < 27) {
        if (str->data[kstr] != cv0[kstr]) {
          exitg11 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg11 = 1;
      }
    } while (exitg11 == 0);
  }

  if (b_bool) {
    kstr = 0;
  } else {
    b_bool = false;
    if (str->size[1] != 28) {
    } else {
      kstr = 0;
      do {
        exitg10 = 0;
        if (kstr + 1 < 29) {
          if (str->data[kstr] != cv1[kstr]) {
            exitg10 = 1;
          } else {
            kstr++;
          }
        } else {
          b_bool = true;
          exitg10 = 1;
        }
      } while (exitg10 == 0);
    }

    if (b_bool) {
      kstr = 1;
    } else {
      b_bool = false;
      if (str->size[1] != 23) {
      } else {
        kstr = 0;
        do {
          exitg9 = 0;
          if (kstr + 1 < 24) {
            if (str->data[kstr] != cv2[kstr]) {
              exitg9 = 1;
            } else {
              kstr++;
            }
          } else {
            b_bool = true;
            exitg9 = 1;
          }
        } while (exitg9 == 0);
      }

      if (b_bool) {
        kstr = 2;
      } else {
        b_bool = false;
        if (str->size[1] != 31) {
        } else {
          kstr = 0;
          do {
            exitg8 = 0;
            if (kstr + 1 < 32) {
              if (str->data[kstr] != cv3[kstr]) {
                exitg8 = 1;
              } else {
                kstr++;
              }
            } else {
              b_bool = true;
              exitg8 = 1;
            }
          } while (exitg8 == 0);
        }

        if (b_bool) {
          kstr = 3;
        } else {
          b_bool = false;
          if (str->size[1] != 28) {
          } else {
            kstr = 0;
            do {
              exitg7 = 0;
              if (kstr + 1 < 29) {
                if (str->data[kstr] != cv4[kstr]) {
                  exitg7 = 1;
                } else {
                  kstr++;
                }
              } else {
                b_bool = true;
                exitg7 = 1;
              }
            } while (exitg7 == 0);
          }

          if (b_bool) {
            kstr = 4;
          } else {
            b_bool = false;
            if (str->size[1] != 29) {
            } else {
              kstr = 0;
              do {
                exitg6 = 0;
                if (kstr + 1 < 30) {
                  if (str->data[kstr] != cv5[kstr]) {
                    exitg6 = 1;
                  } else {
                    kstr++;
                  }
                } else {
                  b_bool = true;
                  exitg6 = 1;
                }
              } while (exitg6 == 0);
            }

            if (b_bool) {
              kstr = 5;
            } else {
              b_bool = false;
              if (str->size[1] != 29) {
              } else {
                kstr = 0;
                do {
                  exitg5 = 0;
                  if (kstr + 1 < 30) {
                    if (str->data[kstr] != cv6[kstr]) {
                      exitg5 = 1;
                    } else {
                      kstr++;
                    }
                  } else {
                    b_bool = true;
                    exitg5 = 1;
                  }
                } while (exitg5 == 0);
              }

              if (b_bool) {
                kstr = 6;
              } else {
                b_bool = false;
                if (str->size[1] != 29) {
                } else {
                  kstr = 0;
                  do {
                    exitg4 = 0;
                    if (kstr + 1 < 30) {
                      if (str->data[kstr] != cv7[kstr]) {
                        exitg4 = 1;
                      } else {
                        kstr++;
                      }
                    } else {
                      b_bool = true;
                      exitg4 = 1;
                    }
                  } while (exitg4 == 0);
                }

                if (b_bool) {
                  kstr = 7;
                } else {
                  b_bool = false;
                  if (str->size[1] != 32) {
                  } else {
                    kstr = 0;
                    do {
                      exitg3 = 0;
                      if (kstr + 1 < 33) {
                        if (str->data[kstr] != cv8[kstr]) {
                          exitg3 = 1;
                        } else {
                          kstr++;
                        }
                      } else {
                        b_bool = true;
                        exitg3 = 1;
                      }
                    } while (exitg3 == 0);
                  }

                  if (b_bool) {
                    kstr = 8;
                  } else {
                    b_bool = false;
                    if (str->size[1] != 30) {
                    } else {
                      kstr = 0;
                      do {
                        exitg2 = 0;
                        if (kstr + 1 < 31) {
                          if (str->data[kstr] != cv9[kstr]) {
                            exitg2 = 1;
                          } else {
                            kstr++;
                          }
                        } else {
                          b_bool = true;
                          exitg2 = 1;
                        }
                      } while (exitg2 == 0);
                    }

                    if (b_bool) {
                      kstr = 9;
                    } else {
                      b_bool = false;
                      if (str->size[1] != 41) {
                      } else {
                        kstr = 0;
                        do {
                          exitg1 = 0;
                          if (kstr + 1 < 42) {
                            if (str->data[kstr] != cv10[kstr]) {
                              exitg1 = 1;
                            } else {
                              kstr++;
                            }
                          } else {
                            b_bool = true;
                            exitg1 = 1;
                          }
                        } while (exitg1 == 0);
                      }

                      if (b_bool) {
                        kstr = 10;
                      } else {
                        kstr = -1;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  switch (kstr) {
   case 0:
    val = (CUSPARSE_POINTER_MODE_HOST);
    break;

   case 1:
    val = (CUSPARSE_POINTER_MODE_DEVICE);
    break;

   case 2:
    val = (CUSPARSE_STATUS_SUCCESS);
    break;

   case 3:
    val = (CUSPARSE_STATUS_NOT_INITIALIZED);
    break;

   case 4:
    val = (CUSPARSE_STATUS_ALLOC_FAILED);
    break;

   case 5:
    val = (CUSPARSE_STATUS_INVALID_VALUE);
    break;

   case 6:
    val = (CUSPARSE_STATUS_ARCH_MISMATCH);
    break;

   case 7:
    val = (CUSPARSE_STATUS_MAPPING_ERROR);
    break;

   case 8:
    val = (CUSPARSE_STATUS_EXECUTION_FAILED);
    break;

   case 9:
    val = (CUSPARSE_STATUS_INTERNAL_ERROR);
    break;

   case 10:
    val = (CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED);
    break;

   default:
    val = -1;
    break;
  }

  return val;
}

void cuSparseGetEnum_initialize(void)
{
}

void cuSparseGetEnum_terminate(void)
{
}

void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxInit_char_T(pEmxArray, numDimensions);
}
