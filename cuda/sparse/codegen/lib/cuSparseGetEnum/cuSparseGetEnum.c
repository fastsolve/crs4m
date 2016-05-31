#include "cuSparseGetEnum.h"
#include "mspack.h"
#include "m2c.h"

int32_T cuSparseGetEnum(const emxArray_char_T *str)
{
  int32_T val;
  boolean_T b_bool;
  int32_T kstr;
  int32_T exitg26;
  static const char_T cv0[26] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'P', 'O', 'I', 'N', 'T', 'E', 'R', '_', 'M', 'O', 'D', 'E', '_', 'H', 'O',
    'S', 'T' };

  int32_T exitg25;
  static const char_T cv1[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'P', 'O', 'I', 'N', 'T', 'E', 'R', '_', 'M', 'O', 'D', 'E', '_', 'D', 'E',
    'V', 'I', 'C', 'E' };

  int32_T exitg24;
  static const char_T cv2[27] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'D', 'I', 'A', 'G', '_', 'T', 'Y', 'P', 'E', '_', 'N', 'O', 'N', '_', 'U',
    'N', 'I', 'T' };

  int32_T exitg23;
  static const char_T cv3[23] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'D', 'I', 'A', 'G', '_', 'T', 'Y', 'P', 'E', '_', 'U', 'N', 'I', 'T' };

  int32_T exitg22;
  static const char_T cv4[24] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'F', 'I', 'L', 'L', '_', 'M', 'O', 'D', 'E', '_', 'L', 'O', 'W', 'E', 'R' };

  int32_T exitg21;
  static const char_T cv5[24] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'F', 'I', 'L', 'L', '_', 'M', 'O', 'D', 'E', '_', 'U', 'P', 'P', 'E', 'R' };

  int32_T exitg20;
  static const char_T cv6[24] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'I', 'N', 'D', 'E', 'X', '_', 'B', 'A', 'S', 'E', '_', 'Z', 'E', 'R', 'O' };

  int32_T exitg19;
  static const char_T cv7[23] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'I', 'N', 'D', 'E', 'X', '_', 'B', 'A', 'S', 'E', '_', 'O', 'N', 'E' };

  int32_T exitg18;
  static const char_T cv8[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'M', 'A', 'T', 'R', 'I', 'X', '_', 'T', 'Y', 'P', 'E', '_', 'G', 'E', 'N',
    'E', 'R', 'A', 'L' };

  int32_T exitg17;
  static const char_T cv9[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'M', 'A', 'T', 'R', 'I', 'X', '_', 'T', 'Y', 'P', 'E', '_', 'S', 'Y', 'M',
    'M', 'E', 'T', 'R', 'I', 'C' };

  int32_T exitg16;
  static const char_T cv10[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'M', 'A', 'T', 'R', 'I', 'X', '_', 'T', 'Y', 'P', 'E', '_', 'H', 'E', 'R',
    'M', 'I', 'T', 'I', 'A', 'N' };

  int32_T exitg15;
  static const char_T cv11[31] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'M', 'A', 'T', 'R', 'I', 'X', '_', 'T', 'Y', 'P', 'E', '_', 'T', 'R', 'I',
    'A', 'N', 'G', 'U', 'L', 'A', 'R' };

  int32_T exitg14;
  static const char_T cv12[22] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'D', 'I', 'R', 'E', 'C', 'T', 'I', 'O', 'N', '_', 'R', 'O', 'W' };

  int32_T exitg13;
  static const char_T cv13[25] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'D', 'I', 'R', 'E', 'C', 'T', 'I', 'O', 'N', '_', 'C', 'O', 'L', 'U', 'M',
    'N' };

  int32_T exitg12;
  static const char_T cv14[32] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'O', 'P', 'E', 'R', 'A', 'T', 'I', 'O', 'N', '_', 'N', 'O', 'N', '_', 'T',
    'R', 'A', 'N', 'S', 'P', 'O', 'S', 'E' };

  int32_T exitg11;
  static const char_T cv15[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'O', 'P', 'E', 'R', 'A', 'T', 'I', 'O', 'N', '_', 'T', 'R', 'A', 'N', 'S',
    'P', 'O', 'S', 'E' };

  int32_T exitg10;
  static const char_T cv16[38] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'O', 'P', 'E', 'R', 'A', 'T', 'I', 'O', 'N', '_', 'C', 'O', 'N', 'J', 'U',
    'G', 'A', 'T', 'E', '_', 'T', 'R', 'A', 'N', 'S', 'P', 'O', 'S', 'E' };

  int32_T exitg9;
  static const char_T cv17[23] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'S', 'U', 'C', 'C', 'E', 'S', 'S' };

  int32_T exitg8;
  static const char_T cv18[31] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'N', 'O', 'T', '_', 'I', 'N', 'I', 'T',
    'I', 'A', 'L', 'I', 'Z', 'E', 'D' };

  int32_T exitg7;
  static const char_T cv19[28] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'A', 'L', 'L', 'O', 'C', '_', 'F', 'A',
    'I', 'L', 'E', 'D' };

  int32_T exitg6;
  static const char_T cv20[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'V', 'A', 'L', 'I', 'D', '_',
    'V', 'A', 'L', 'U', 'E' };

  int32_T exitg5;
  static const char_T cv21[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'A', 'R', 'C', 'H', '_', 'M', 'I', 'S',
    'M', 'A', 'T', 'C', 'H' };

  int32_T exitg4;
  static const char_T cv22[29] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'P', 'P', 'I', 'N', 'G', '_',
    'E', 'R', 'R', 'O', 'R' };

  int32_T exitg3;
  static const char_T cv23[32] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'E', 'X', 'E', 'C', 'U', 'T', 'I', 'O',
    'N', '_', 'F', 'A', 'I', 'L', 'E', 'D' };

  int32_T exitg2;
  static const char_T cv24[30] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'I', 'N', 'T', 'E', 'R', 'N', 'A', 'L',
    '_', 'E', 'R', 'R', 'O', 'R' };

  int32_T exitg1;
  static const char_T cv25[41] = { 'C', 'U', 'S', 'P', 'A', 'R', 'S', 'E', '_',
    'S', 'T', 'A', 'T', 'U', 'S', '_', 'M', 'A', 'T', 'R', 'I', 'X', '_', 'T',
    'Y', 'P', 'E', '_', 'N', 'O', 'T', '_', 'S', 'U', 'P', 'P', 'O', 'R', 'T',
    'E', 'D' };

  b_bool = false;
  if (str->size[1] != 26) {
  } else {
    kstr = 0;
    do {
      exitg26 = 0;
      if (kstr + 1 < 27) {
        if (str->data[kstr] != cv0[kstr]) {
          exitg26 = 1;
        } else {
          kstr++;
        }
      } else {
        b_bool = true;
        exitg26 = 1;
      }
    } while (exitg26 == 0);
  }

  if (b_bool) {
    kstr = 0;
  } else {
    b_bool = false;
    if (str->size[1] != 28) {
    } else {
      kstr = 0;
      do {
        exitg25 = 0;
        if (kstr + 1 < 29) {
          if (str->data[kstr] != cv1[kstr]) {
            exitg25 = 1;
          } else {
            kstr++;
          }
        } else {
          b_bool = true;
          exitg25 = 1;
        }
      } while (exitg25 == 0);
    }

    if (b_bool) {
      kstr = 1;
    } else {
      b_bool = false;
      if (str->size[1] != 27) {
      } else {
        kstr = 0;
        do {
          exitg24 = 0;
          if (kstr + 1 < 28) {
            if (str->data[kstr] != cv2[kstr]) {
              exitg24 = 1;
            } else {
              kstr++;
            }
          } else {
            b_bool = true;
            exitg24 = 1;
          }
        } while (exitg24 == 0);
      }

      if (b_bool) {
        kstr = 2;
      } else {
        b_bool = false;
        if (str->size[1] != 23) {
        } else {
          kstr = 0;
          do {
            exitg23 = 0;
            if (kstr + 1 < 24) {
              if (str->data[kstr] != cv3[kstr]) {
                exitg23 = 1;
              } else {
                kstr++;
              }
            } else {
              b_bool = true;
              exitg23 = 1;
            }
          } while (exitg23 == 0);
        }

        if (b_bool) {
          kstr = 3;
        } else {
          b_bool = false;
          if (str->size[1] != 24) {
          } else {
            kstr = 0;
            do {
              exitg22 = 0;
              if (kstr + 1 < 25) {
                if (str->data[kstr] != cv4[kstr]) {
                  exitg22 = 1;
                } else {
                  kstr++;
                }
              } else {
                b_bool = true;
                exitg22 = 1;
              }
            } while (exitg22 == 0);
          }

          if (b_bool) {
            kstr = 4;
          } else {
            b_bool = false;
            if (str->size[1] != 24) {
            } else {
              kstr = 0;
              do {
                exitg21 = 0;
                if (kstr + 1 < 25) {
                  if (str->data[kstr] != cv5[kstr]) {
                    exitg21 = 1;
                  } else {
                    kstr++;
                  }
                } else {
                  b_bool = true;
                  exitg21 = 1;
                }
              } while (exitg21 == 0);
            }

            if (b_bool) {
              kstr = 5;
            } else {
              b_bool = false;
              if (str->size[1] != 24) {
              } else {
                kstr = 0;
                do {
                  exitg20 = 0;
                  if (kstr + 1 < 25) {
                    if (str->data[kstr] != cv6[kstr]) {
                      exitg20 = 1;
                    } else {
                      kstr++;
                    }
                  } else {
                    b_bool = true;
                    exitg20 = 1;
                  }
                } while (exitg20 == 0);
              }

              if (b_bool) {
                kstr = 6;
              } else {
                b_bool = false;
                if (str->size[1] != 23) {
                } else {
                  kstr = 0;
                  do {
                    exitg19 = 0;
                    if (kstr + 1 < 24) {
                      if (str->data[kstr] != cv7[kstr]) {
                        exitg19 = 1;
                      } else {
                        kstr++;
                      }
                    } else {
                      b_bool = true;
                      exitg19 = 1;
                    }
                  } while (exitg19 == 0);
                }

                if (b_bool) {
                  kstr = 7;
                } else {
                  b_bool = false;
                  if (str->size[1] != 28) {
                  } else {
                    kstr = 0;
                    do {
                      exitg18 = 0;
                      if (kstr + 1 < 29) {
                        if (str->data[kstr] != cv8[kstr]) {
                          exitg18 = 1;
                        } else {
                          kstr++;
                        }
                      } else {
                        b_bool = true;
                        exitg18 = 1;
                      }
                    } while (exitg18 == 0);
                  }

                  if (b_bool) {
                    kstr = 8;
                  } else {
                    b_bool = false;
                    if (str->size[1] != 30) {
                    } else {
                      kstr = 0;
                      do {
                        exitg17 = 0;
                        if (kstr + 1 < 31) {
                          if (str->data[kstr] != cv9[kstr]) {
                            exitg17 = 1;
                          } else {
                            kstr++;
                          }
                        } else {
                          b_bool = true;
                          exitg17 = 1;
                        }
                      } while (exitg17 == 0);
                    }

                    if (b_bool) {
                      kstr = 9;
                    } else {
                      b_bool = false;
                      if (str->size[1] != 30) {
                      } else {
                        kstr = 0;
                        do {
                          exitg16 = 0;
                          if (kstr + 1 < 31) {
                            if (str->data[kstr] != cv10[kstr]) {
                              exitg16 = 1;
                            } else {
                              kstr++;
                            }
                          } else {
                            b_bool = true;
                            exitg16 = 1;
                          }
                        } while (exitg16 == 0);
                      }

                      if (b_bool) {
                        kstr = 10;
                      } else {
                        b_bool = false;
                        if (str->size[1] != 31) {
                        } else {
                          kstr = 0;
                          do {
                            exitg15 = 0;
                            if (kstr + 1 < 32) {
                              if (str->data[kstr] != cv11[kstr]) {
                                exitg15 = 1;
                              } else {
                                kstr++;
                              }
                            } else {
                              b_bool = true;
                              exitg15 = 1;
                            }
                          } while (exitg15 == 0);
                        }

                        if (b_bool) {
                          kstr = 11;
                        } else {
                          b_bool = false;
                          if (str->size[1] != 22) {
                          } else {
                            kstr = 0;
                            do {
                              exitg14 = 0;
                              if (kstr + 1 < 23) {
                                if (str->data[kstr] != cv12[kstr]) {
                                  exitg14 = 1;
                                } else {
                                  kstr++;
                                }
                              } else {
                                b_bool = true;
                                exitg14 = 1;
                              }
                            } while (exitg14 == 0);
                          }

                          if (b_bool) {
                            kstr = 12;
                          } else {
                            b_bool = false;
                            if (str->size[1] != 25) {
                            } else {
                              kstr = 0;
                              do {
                                exitg13 = 0;
                                if (kstr + 1 < 26) {
                                  if (str->data[kstr] != cv13[kstr]) {
                                    exitg13 = 1;
                                  } else {
                                    kstr++;
                                  }
                                } else {
                                  b_bool = true;
                                  exitg13 = 1;
                                }
                              } while (exitg13 == 0);
                            }

                            if (b_bool) {
                              kstr = 13;
                            } else {
                              b_bool = false;
                              if (str->size[1] != 32) {
                              } else {
                                kstr = 0;
                                do {
                                  exitg12 = 0;
                                  if (kstr + 1 < 33) {
                                    if (str->data[kstr] != cv14[kstr]) {
                                      exitg12 = 1;
                                    } else {
                                      kstr++;
                                    }
                                  } else {
                                    b_bool = true;
                                    exitg12 = 1;
                                  }
                                } while (exitg12 == 0);
                              }

                              if (b_bool) {
                                kstr = 14;
                              } else {
                                b_bool = false;
                                if (str->size[1] != 28) {
                                } else {
                                  kstr = 0;
                                  do {
                                    exitg11 = 0;
                                    if (kstr + 1 < 29) {
                                      if (str->data[kstr] != cv15[kstr]) {
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
                                  kstr = 15;
                                } else {
                                  b_bool = false;
                                  if (str->size[1] != 38) {
                                  } else {
                                    kstr = 0;
                                    do {
                                      exitg10 = 0;
                                      if (kstr + 1 < 39) {
                                        if (str->data[kstr] != cv16[kstr]) {
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
                                    kstr = 16;
                                  } else {
                                    b_bool = false;
                                    if (str->size[1] != 23) {
                                    } else {
                                      kstr = 0;
                                      do {
                                        exitg9 = 0;
                                        if (kstr + 1 < 24) {
                                          if (str->data[kstr] != cv17[kstr]) {
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
                                      kstr = 17;
                                    } else {
                                      b_bool = false;
                                      if (str->size[1] != 31) {
                                      } else {
                                        kstr = 0;
                                        do {
                                          exitg8 = 0;
                                          if (kstr + 1 < 32) {
                                            if (str->data[kstr] != cv18[kstr]) {
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
                                        kstr = 18;
                                      } else {
                                        b_bool = false;
                                        if (str->size[1] != 28) {
                                        } else {
                                          kstr = 0;
                                          do {
                                            exitg7 = 0;
                                            if (kstr + 1 < 29) {
                                              if (str->data[kstr] != cv19[kstr])
                                              {
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
                                          kstr = 19;
                                        } else {
                                          b_bool = false;
                                          if (str->size[1] != 29) {
                                          } else {
                                            kstr = 0;
                                            do {
                                              exitg6 = 0;
                                              if (kstr + 1 < 30) {
                                                if (str->data[kstr] != cv20[kstr])
                                                {
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
                                            kstr = 20;
                                          } else {
                                            b_bool = false;
                                            if (str->size[1] != 29) {
                                            } else {
                                              kstr = 0;
                                              do {
                                                exitg5 = 0;
                                                if (kstr + 1 < 30) {
                                                  if (str->data[kstr] !=
                                                      cv21[kstr]) {
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
                                              kstr = 21;
                                            } else {
                                              b_bool = false;
                                              if (str->size[1] != 29) {
                                              } else {
                                                kstr = 0;
                                                do {
                                                  exitg4 = 0;
                                                  if (kstr + 1 < 30) {
                                                    if (str->data[kstr] !=
                                                        cv22[kstr]) {
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
                                                kstr = 22;
                                              } else {
                                                b_bool = false;
                                                if (str->size[1] != 32) {
                                                } else {
                                                  kstr = 0;
                                                  do {
                                                    exitg3 = 0;
                                                    if (kstr + 1 < 33) {
                                                      if (str->data[kstr] !=
                                                          cv23[kstr]) {
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
                                                  kstr = 23;
                                                } else {
                                                  b_bool = false;
                                                  if (str->size[1] != 30) {
                                                  } else {
                                                    kstr = 0;
                                                    do {
                                                      exitg2 = 0;
                                                      if (kstr + 1 < 31) {
                                                        if (str->data[kstr] !=
                                                            cv24[kstr]) {
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
                                                    kstr = 24;
                                                  } else {
                                                    b_bool = false;
                                                    if (str->size[1] != 41) {
                                                    } else {
                                                      kstr = 0;
                                                      do {
                                                        exitg1 = 0;
                                                        if (kstr + 1 < 42) {
                                                          if (str->data[kstr] !=
                                                              cv25[kstr]) {
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
                                                      kstr = 25;
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
    val = (CUSPARSE_DIAG_TYPE_NON_UNIT);
    break;

   case 3:
    val = (CUSPARSE_DIAG_TYPE_UNIT);
    break;

   case 4:
    val = (CUSPARSE_FILL_MODE_LOWER);
    break;

   case 5:
    val = (CUSPARSE_FILL_MODE_UPPER);
    break;

   case 6:
    val = (CUSPARSE_INDEX_BASE_ZERO);
    break;

   case 7:
    val = (CUSPARSE_INDEX_BASE_ONE);
    break;

   case 8:
    val = (CUSPARSE_MATRIX_TYPE_GENERAL);
    break;

   case 9:
    val = (CUSPARSE_MATRIX_TYPE_SYMMETRIC);
    break;

   case 10:
    val = (CUSPARSE_MATRIX_TYPE_HERMITIAN);
    break;

   case 11:
    val = (CUSPARSE_MATRIX_TYPE_TRIANGULAR);
    break;

   case 12:
    val = (CUSPARSE_DIRECTION_ROW);
    break;

   case 13:
    val = (CUSPARSE_DIRECTION_COLUMN);
    break;

   case 14:
    val = (CUSPARSE_OPERATION_NON_TRANSPOSE);
    break;

   case 15:
    val = (CUSPARSE_OPERATION_TRANSPOSE);
    break;

   case 16:
    val = (CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE);
    break;

   case 17:
    val = (CUSPARSE_STATUS_SUCCESS);
    break;

   case 18:
    val = (CUSPARSE_STATUS_NOT_INITIALIZED);
    break;

   case 19:
    val = (CUSPARSE_STATUS_ALLOC_FAILED);
    break;

   case 20:
    val = (CUSPARSE_STATUS_INVALID_VALUE);
    break;

   case 21:
    val = (CUSPARSE_STATUS_ARCH_MISMATCH);
    break;

   case 22:
    val = (CUSPARSE_STATUS_MAPPING_ERROR);
    break;

   case 23:
    val = (CUSPARSE_STATUS_EXECUTION_FAILED);
    break;

   case 24:
    val = (CUSPARSE_STATUS_INTERNAL_ERROR);
    break;

   case 25:
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

void emxInitArray_char_T(emxArray_char_T **pEmxArray, int32_T numDimensions)
{
  emxInit_char_T(pEmxArray, numDimensions);
}
