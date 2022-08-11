#include "fft_processor_spqlios.h"

#include <cstdlib>
#include <cassert>

namespace {
  typedef struct {
    uint64_t n;
    double *aligned_trig_tables;
    double *aligned_data;
    void *buf;
  } FFT_PRECOMP;

  typedef struct {
    uint64_t n;
    double *aligned_trig_tables;
    double *aligned_data;
    void *buf;
  } IFFT_PRECOMP;

  //trig_tables:
  //|cos0,cos1,cos2,cos3|sin0,sin1,sin2,sin3|cos4,cos5,cos6,cos7|.... -> n/4
  //|cos0,cos2,cos4,cos6|sin0,sin2,sin4,sin6|cos8,cos10,cos12,cos14|.... -> n/4
  //...
  //|cos0,cos2,cos4,cos6|sin0,sin2,sin4,sin6|cos8,cos10,cos12,cos14|.... -> n/4


  //trig_tables:
  //|cos0,cos1,cos2,cos3|sin0,sin1,sin2,sin3|cos4,cos5,cos6,cos7|.... -> n/4
  //|cos0,cos2,cos4,cos6|sin0,sin2,sin4,sin6|cos8,cos10,cos12,cos14|.... -> n/4
  //...
  //|cos0,cos2,cos4,cos6|sin0,sin2,sin4,sin6|cos8,cos10,cos12,cos14|.... -> n/4

  inline void dotp4(double * res, const double * a, const double * b) {
    for (int32_t i = 0; i < 4; i++) res[i] = a[i] * b[i];
  }

  inline void add4(double * res, const double * a, const double * b) {
    for (int32_t i = 0; i < 4; i++) res[i] = a[i] + b[i];
  }

  inline void sub4(double * res, const double * a, const double * b) {
    for (int32_t i = 0; i < 4; i++) res[i] = a[i] - b[i];
  }

  inline void copy4(double * res, const double * a) {
    for (int32_t i = 0; i < 4; i++) res[i] = a[i];
  }


  double accurate_cos(int32_t i, int32_t n) { //cos(2pi*i/n)
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4) return cos(2. * M_PI * (n - i) / ((double) n));
    if (i >= 2 * n / 4) return -cos(2. * M_PI * (i - n / 2) / ((double) n));
    if (i >= 1 * n / 4) return -cos(2. * M_PI * (n / 2 - i) / ((double) n));
    return cos(2. * M_PI * (i) / ((double) n));
  }

  double accurate_sin(int32_t i, int32_t n) { //sin(2pi*i/n)
    i = ((i % n) + n) % n;
    if (i >= 3 * n / 4) return -sin(2. * M_PI * (n - i) / ((double) n));
    if (i >= 2 * n / 4) return -sin(2. * M_PI * (i - n / 2) / ((double) n));
    if (i >= 1 * n / 4) return sin(2. * M_PI * (n / 2 - i) / ((double) n));
    return sin(2. * M_PI * (i) / ((double) n));
  }
}

extern "C" void *new_fft_table(int32_t nn) {
  assert(nn >= 16);
  assert((nn & (nn - 1)) == 0);
  int32_t n = 2 * nn;
  int32_t ns4 = n / 4;
  FFT_PRECOMP *reps = (FFT_PRECOMP *) std::aligned_alloc(64,sizeof(FFT_PRECOMP));
  void *buf = std::aligned_alloc(64,n * 8 + nn * 8);
  // uint64_t aligned_addr = (((uint64_t) buf) + 31) & 0xFFFFFFFFFFFFFFE0l;
  uint64_t aligned_addr = (uint64_t) buf;
  reps->n = n;
  reps->aligned_trig_tables = (double *) aligned_addr;
  reps->aligned_data = (double *) (aligned_addr + n * 8);
  reps->buf = buf;
  double *ptr = reps->aligned_trig_tables;
  // first iteration
  const int32_t j_0 = n / 8;
  for (int32_t k = 0; k < 4; k++)
    *(ptr++) = accurate_cos(-j_0 * k, n);
  for (int32_t k = 0; k < 4; k++)
    *(ptr++) = accurate_sin(-j_0 * k, n);
  // subsequent iterations
  for (int32_t halfnn = 8; halfnn < ns4; halfnn *= 2) {
    int32_t nn = 2 * halfnn;
    int32_t j = n / nn;
    //cerr << "- b: " << halfnn  << "(offset: " << (ptr-reps->trig_tables) << ", mult: " << j << ")" << endl;
    for (int32_t i = 0; i < halfnn; i += 8) {
      //cerr << "--- i: " << i << endl;
      for (int32_t k = 0; k < 8; k++)
        *(ptr++) = accurate_cos(-j * (i + k), n);
      for (int32_t k = 0; k < 8; k++)
        *(ptr++) = accurate_sin(-j * (i + k), n);
    }
  }
  //last iteration
  for (int32_t i = 0; i < ns4; i += 8) {
    for (int32_t k = 0; k < 8; k++)
      *(ptr++) = accurate_cos(-(i + k), n);
    for (int32_t k = 0; k < 8; k++)
      *(ptr++) = accurate_sin(-(i + k), n);
  }
  return reps;
}

extern "C" double *fft_table_get_buffer(const void *tables) {
  FFT_PRECOMP *reps = (FFT_PRECOMP *) tables;
  return reps->aligned_data;
}

extern "C" double *ifft_table_get_buffer(const void *tables) {
  IFFT_PRECOMP *reps = (IFFT_PRECOMP *) tables;
  return reps->aligned_data;
}

//c has size n/2
extern "C" void fft_model(const void *tables) {
  double tmp0[4];
  double tmp1[4];
  double tmp2[4];
  double tmp3[4];
  FFT_PRECOMP *fft_tables = (FFT_PRECOMP *) tables;
  const int32_t n = fft_tables->n;
  const double *trig_tables = fft_tables->aligned_trig_tables;
  double *c = fft_tables->aligned_data;

  int32_t ns4 = n / 4;
  double *pre = c;   //size n/4
  double *pim = c + ns4; //size n/4

  //general loop
  //size 2
  {
    //[1  1]
    //[1 -1]
    //   [1  1]
    //   [1 -1]
    for (int32_t block = 0; block < ns4; block += 4) {
      double *d0 = pre + block;
      double *d1 = pim + block;
      tmp0[0] = d0[0];
      tmp0[1] = d0[0];
      tmp0[2] = d0[2];
      tmp0[3] = d0[2];
      tmp1[0] = d0[1];
      tmp1[1] = -d0[1];
      tmp1[2] = d0[3];
      tmp1[3] = -d0[3];
      add4(d0, tmp0, tmp1);
      tmp0[0] = d1[0];
      tmp0[1] = d1[0];
      tmp0[2] = d1[2];
      tmp0[3] = d1[2];
      tmp1[0] = d1[1];
      tmp1[1] = -d1[1];
      tmp1[2] = d1[3];
      tmp1[3] = -d1[3];
      add4(d1, tmp0, tmp1);
    }
  }

  //size 4
  //[1  0  1  0]
  //[0  1  0 -i]
  //[1  0 -1  0]
  //[0  1  0  i]
  // r0 + r2  i0 + i2
  // r1 + i3  i1 - r3
  // r0 - r2  i0 - i2
  // r1 - i3  i1 + r3
  {
    for (int32_t block = 0; block < ns4; block += 4) {
      double *re = pre + block;
      double *im = pim + block;
      tmp0[0] = re[0];
      tmp0[1] = re[1];
      tmp0[2] = re[0];
      tmp0[3] = re[1];
      tmp1[0] = re[2];
      tmp1[1] = im[3];
      tmp1[2] = -re[2];
      tmp1[3] = -im[3];
      tmp2[0] = im[0];
      tmp2[1] = im[1];
      tmp2[2] = im[0];
      tmp2[3] = im[1];
      tmp3[0] = im[2];
      tmp3[1] = -re[3];
      tmp3[2] = -im[2];
      tmp3[3] = re[3];
      add4(re, tmp0, tmp1);
      add4(im, tmp2, tmp3);
    }
  }

  //general loop
  const double *cur_tt = trig_tables;
  for (int32_t halfnn = 4; halfnn < ns4; halfnn *= 2) {
    int32_t nn = 2 * halfnn;
    for (int32_t block = 0; block < ns4; block += nn) {
      for (int32_t off = 0; off < halfnn; off += 4) {
        double *re0 = pre + block + off;
        double *im0 = pim + block + off;
        double *re1 = pre + block + halfnn + off;
        double *im1 = pim + block + halfnn + off;
        const double *tcs = cur_tt + 2 * off;
        const double *tsn = tcs + 4;
        dotp4(tmp0, re1, tcs); // re*cos
        dotp4(tmp1, re1, tsn); // re*sin
        dotp4(tmp2, im1, tcs); // im*cos
        dotp4(tmp3, im1, tsn); // im*sin
        sub4(tmp0, tmp0, tmp3); // re2
        add4(tmp1, tmp1, tmp2); // im2
        add4(tmp2, re0, tmp0); // re + re
        add4(tmp3, im0, tmp1); // im + im
        sub4(tmp0, re0, tmp0); // re - re
        sub4(tmp1, im0, tmp1); // im - im
        copy4(re0, tmp2);
        copy4(im0, tmp3);
        copy4(re1, tmp0);
        copy4(im1, tmp1);
      }
    }
    cur_tt += nn;
  }

  //multiply by omb^j
  for (int32_t j = 0; j < ns4; j += 4) {
    const double *r0 = cur_tt + 2 * j;
    const double *r1 = r0 + 4;
    //(re*cos-im*sin) + i (im*cos+re*sin)
    double *d0 = pre + j;
    double *d1 = pim + j;
    dotp4(tmp0, d0, r0); //re*cos
    dotp4(tmp1, d1, r0); //im*cos
    dotp4(tmp2, d0, r1); //re*sin
    dotp4(tmp3, d1, r1); //im*sin
    sub4(d0, tmp0, tmp3);
    add4(d1, tmp1, tmp2);
  }
}

extern "C" void *new_ifft_table(int32_t nn) {
  assert(nn >= 16);
  assert((nn & (nn - 1)) == 0);
  int32_t n = 2 * nn;
  int32_t ns4 = n / 4;
  IFFT_PRECOMP *reps = (IFFT_PRECOMP *) std::aligned_alloc(64,sizeof(IFFT_PRECOMP));
  void *buf = std::aligned_alloc(64,n * 8 + nn * 8);
  // uint64_t aligned_addr = (((uint64_t) buf) + 31) & 0xFFFFFFFFFFFFFFE0l;
  uint64_t aligned_addr = (uint64_t) buf;
  reps->n = n;
  reps->aligned_trig_tables = (double *) aligned_addr;
  reps->aligned_data = (double *) (aligned_addr + n * 8);
  reps->buf = buf;
  reps->n = n;
  double *ptr = reps->aligned_trig_tables;
  //first iteration
  for (int32_t j = 0; j < ns4; j += 8) {
    for (int32_t k = 0; k < 8; k++)
      *(ptr++) = accurate_cos(j + k, n);
    for (int32_t k = 0; k < 8; k++)
      *(ptr++) = accurate_sin(j + k, n);
  }
  //subsequent iterations
  for (int32_t nn = ns4; nn >= 16; nn /= 2) {
    int32_t halfnn = nn / 2;
    int32_t j = n / nn;
    for (int32_t i = 0; i < halfnn; i += 8) {
      for (int32_t k = 0; k < 8; k++)
        *(ptr++) = accurate_cos(j * (i + k), n);
      for (int32_t k = 0; k < 8; k++)
        *(ptr++) = accurate_sin(j * (i + k), n);
    }
  }
  // last iteration
  const int32_t j_last = n / 8;
  for (int32_t k = 0; k < 4; k++)
    *(ptr++) = accurate_cos(j_last * k, n);
  for (int32_t k = 0; k < 4; k++)
    *(ptr++) = accurate_sin(j_last * k, n);
  return reps;
}

//c has size n/2
extern "C" void ifft_model(void *tables) {
  double tmp0[4];
  double tmp1[4];
  double tmp2[4];
  double tmp3[4];
  IFFT_PRECOMP *fft_tables = (IFFT_PRECOMP *) tables;
  const int32_t n = fft_tables->n;
  const double *trig_tables = fft_tables->aligned_trig_tables;
  double *c = fft_tables->aligned_data;

  int32_t ns4 = n / 4;
  double *are = c;  //size n/4
  double *aim = c + ns4; //size n/4

  //multiply by omega^j
  for (int32_t j = 0; j < ns4; j += 4) {
    const double *r0 = trig_tables + 2 * j;
    const double *r1 = r0 + 4;
    //(re*cos-im*sin) + i (im*cos+re*sin)
    double *d0 = are + j;
    double *d1 = aim + j;
    dotp4(tmp0, d0, r0); //re*cos
    dotp4(tmp1, d1, r0); //im*cos
    dotp4(tmp2, d0, r1); //re*sin
    dotp4(tmp3, d1, r1); //im*sin
    sub4(d0, tmp0, tmp3);
    add4(d1, tmp1, tmp2);
  }


  //at the beginning of iteration nn
  // a_{j,i} has P_{i%nn}(omega^j) 
  // where j between [rev(1) and rev(3)[
  // and i between [0 and nn[
  const double *cur_tt = trig_tables;
  for (int32_t nn = ns4; nn >= 8; nn /= 2) {
    int32_t halfnn = nn / 2;
    cur_tt += 2 * nn;
    for (int32_t block = 0; block < ns4; block += nn) {
      for (int32_t off = 0; off < halfnn; off += 4) {
        double *d00 = are + block + off;
        double *d01 = aim + block + off;
        double *d10 = are + block + halfnn + off;
        double *d11 = aim + block + halfnn + off;
        add4(tmp0, d00, d10); // re + re
        add4(tmp1, d01, d11); // im + im
        sub4(tmp2, d00, d10); // re - re
        sub4(tmp3, d01, d11); // im - im
        copy4(d00, tmp0);
        copy4(d01, tmp1);
        const double *r0 = cur_tt + 2 * off;
        const double *r1 = r0 + 4;
        dotp4(tmp0, tmp2, r0); //re*cos
        dotp4(tmp1, tmp3, r1); //im*sin
        sub4(d10, tmp0, tmp1);
        dotp4(tmp0, tmp2, r1); //re*sin
        dotp4(tmp1, tmp3, r0); //im*cos
        add4(d11, tmp0, tmp1);
      }
    }
  }

  //size 4
  {
    for (int32_t block = 0; block < ns4; block += 4) {
      double *d0 = are + block;
      double *d1 = aim + block;
      tmp0[0] = d0[0];
      tmp0[1] = d0[1];
      tmp0[2] = d0[0];
      tmp0[3] = -d1[1];
      tmp1[0] = d0[2];
      tmp1[1] = d0[3];
      tmp1[2] = -d0[2];
      tmp1[3] = d1[3];
      tmp2[0] = d1[0];
      tmp2[1] = d1[1];
      tmp2[2] = d1[0];
      tmp2[3] = d0[1];
      tmp3[0] = d1[2];
      tmp3[1] = d1[3];
      tmp3[2] = -d1[2];
      tmp3[3] = -d0[3];
      add4(d0, tmp0, tmp1);
      add4(d1, tmp2, tmp3);
    }
  }

  //size 2
  {
    for (int32_t block = 0; block < ns4; block += 4) {
      double *d0 = are + block;
      double *d1 = aim + block;
      tmp0[0] = d0[0];
      tmp0[1] = d0[0];
      tmp0[2] = d0[2];
      tmp0[3] = d0[2];
      tmp1[0] = d0[1];
      tmp1[1] = -d0[1];
      tmp1[2] = d0[3];
      tmp1[3] = -d0[3];
      add4(d0, tmp0, tmp1);
      tmp0[0] = d1[0];
      tmp0[1] = d1[0];
      tmp0[2] = d1[2];
      tmp0[3] = d1[2];
      tmp1[0] = d1[1];
      tmp1[1] = -d1[1];
      tmp1[2] = d1[3];
      tmp1[3] = -d1[3];
      add4(d1, tmp0, tmp1);
    }
  }
}

