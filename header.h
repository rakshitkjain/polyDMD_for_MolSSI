//MD: header.h C++ and library headers (Revision Date: 31st March 2021)
//Global headers
#ifndef _HEADER_H
#define _HEADER_H
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <random>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <ctime>
#define PI 3.14159265
#ifdef _OPENMP
  #include <omp.h>
  #define _IS_SERIAL_ 0
#else
  #define omp_get_num_threads() 1
  #define omp_get_thread_num() 0
  #define _IS_SERIAL_ 1
#endif
using namespace std;
#endif
