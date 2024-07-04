#include "fft_processor_mkl.hpp"

// FFT_Processor_MKL is thread-safe
thread_local FFT_Processor_MKL<TFHEpp::lvl1param::n> fftplvl1;
thread_local FFT_Processor_MKL<TFHEpp::lvl2param::n> fftplvl2;