// This file contains an implementation of a vectorized cosine, which
// is based in part on the implementations in the library "SLEEF" by
// Naoki Shibata. SLEEF was used under the Boost Software License,
// Version 1.0. The original source file contained the following
// copyright notice:
//
//   //          Copyright Naoki Shibata 2010 - 2018.
//   // Distributed under the Boost Software License, Version 1.0.
//   //    (See accompanying file LICENSE.txt or copy at
//   //          http://www.boost.org/LICENSE_1_0.txt)
//
// SLEEF was used under the following license, which is not necessarily the license that applies to this file:
//
//         Boost Software License - Version 1.0 - August 17th, 2003
//         
//         Permission is hereby granted, free of charge, to any person or organization
//         obtaining a copy of the software and accompanying documentation covered by
//         this license (the "Software") to use, reproduce, display, distribute,
//         execute, and transmit the Software, and to prepare derivative works of the
//         Software, and to permit third-parties to whom the Software is furnished to
//         do so, all subject to the following:
//         
//         The copyright notices in the Software and this entire statement, including
//         the above license grant, this restriction and the following disclaimer,
//         must be included in all copies of the Software, in whole or in part, and
//         all derivative works of the Software, unless such copies or derivative
//         works are solely in the form of machine-executable object code generated by
//         a source language processor.
//         
//         THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//         IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//         FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
//         SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
//         FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
//         ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//         DEALINGS IN THE SOFTWARE.


#include "crpropa/magneticField/TD13Field.h"
#include "crpropa/Units.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "kiss/logger.h"

#include <iostream>

#ifdef FAST_TD13
#include <immintrin.h>
#endif

namespace crpropa {

  std::vector<double> logspace(double start, double stop, size_t N) {

    double delta = stop - start;
    std::vector<double> values = std::vector<double>(N, 0.);
    for (int i=0; i<N; i++) {
      values[i] = pow(10, ((double) i) / ((double) (N-1)) * delta + start);
    }
    return values;
  }

#ifdef FAST_TD13

  // this is the second part of hsum_double_avx
  double hsum_double_sse(__m128d v) {
    __m128d high64 = _mm_unpackhi_pd(v, v);
    return _mm_cvtsd_f64(_mm_add_sd(v, high64));
  }
#endif // defined(FAST_TD13)

  TD13Field::TD13Field(double Brms, double kmin, double kmax, double s, int Nm, int seed) {

    // NOTE: the use of the turbulence bend-over scale in the TD13 paper is quite confusing to
    // me. The paper states that k = l_0 * <k tilde> would be used throughout, yet
    // surely they do not mean to say that l_0 * <k tilde> should be used for the k in the
    // scalar product in eq. 2? In this implementation, I've only multiplied in the l_0
    // in the computation of the Gk, not the actual <k>s used for planar wave evaluation,
    // since this would yield obviously wrong results...

#ifdef FAST_TD13
    KISS_LOG_INFO << "TD13Field: Using SIMD TD13 implementation" << std::endl;

    // In principle, we could dynamically dispatch to the non-SIMD version in
    // this case. However, this complicates the code, incurs runtime overhead,
    // and is unlikely to happen since SSE3 is quite well supported.
    // TODO: this is currently uncommented b/c sleef seems to fail to provide
    // the cpuid function
    //if (!check_sse()) {
    //  throw std::runtime_error("TD13Field: This code was compiled with SIMD support (SSE1-3), but it is running on a CPU that does not support these instructions. Please set USE_SIMD to OFF in CMake and recompile CRPropa.");
    //}
#endif

    if (kmin > kmax) {
      throw std::runtime_error("TD13Field: kmin > kmax");
    }

    if (Nm <= 1) {
      throw std::runtime_error("TD13Field: Nm <= 1. We need at least two wavemodes in order to generate the k distribution properly, and besides -- *what are you doing?!*");
    }

    if (kmin <= 0) {
      throw std::runtime_error("TD13Field: kmin <= 0");
    } 

    Random random;
    if (seed != 0) { // copied from initTurbulence
      random.seed(seed);
    }

    // initialize everything
    this->Nm = Nm;

    xi = std::vector<Vector3d>(Nm, Vector3d(0.));
    kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
    phi = std::vector<double>(Nm, 0.);
    costheta = std::vector<double>(Nm, 0.);
    beta = std::vector<double>(Nm, 0.);
    Ak = std::vector<double>(Nm, 0.);

    k = logspace(log10(kmin), log10(kmax), Nm);

    // compute Ak
    double delta_k0 = (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]
    //on second thought, this is probably unnecessary since it's just a factor and will get
    //normalized out anyways.

    double Ak2_sum = 0; // sum of Ak^2 over all k
    //for this loop, the Ak array actually contains Gk*delta_k (ie non-normalized Ak^2)
    for (int i=0; i<Nm; i++) {
      double k = this->k[i];
      double Gk = pow(k, -s);
      Ak[i] = Gk * delta_k0 * k;
      Ak2_sum += Ak[i];
    }
    //only in this loop are the actual Ak computed and stored
    //(this two-step process is necessary in order to normalize the values properly)
    for (int i=0; i<Nm; i++) {
      Ak[i] = sqrt(Ak[i] / Ak2_sum * 2) * Brms;
    }

    // generate direction, phase, and polarization for each wavemode
    for (int i=0; i<Nm; i++) {
      // phi, costheta, and sintheta are for drawing vectors with
      // uniform distribution on the unit sphere.
      // This is similar to Random::randVector(): their t is our phi,
      // z is costheta, and r is sintheta. Our kappa is equivalent to
      // the return value of randVector(); however, TD13 then reuse
      // these values to generate a random vector perpendicular to kappa.
      double phi = random.randUniform(-M_PI, M_PI);
      double costheta = random.randUniform(-1., 1.);
      double sintheta = sqrt(1 - costheta*costheta);

      double alpha = random.randUniform(0, 2*M_PI);
      double beta = random.randUniform(0, 2*M_PI);

      Vector3d kappa = Vector3d ( sintheta * cos(phi), sintheta*sin(phi), costheta );
      Vector3d xi = Vector3d ( costheta*cos(phi)*cos(alpha) + sin(phi)*sin(alpha),
			       costheta*sin(phi)*cos(alpha) - cos(phi)*sin(alpha),
			       -sintheta*cos(alpha) );

      this->xi[i] = xi;
      this->kappa[i] = kappa;
      this->phi[i] = phi;
      this->costheta[i] = costheta;
      this->beta[i] = beta;
    }
    //copy data into AVX-compatible arrays
    avx_Nm = ( (Nm + 4 - 1)/4 ) * 4; //round up to next larger multiple of 4: align is 256 = 4 * sizeof(double) bit
    avx_data = std::vector<double>(itotal*avx_Nm + 3, 0.);

    //get the first 256-bit aligned element
    size_t size = avx_data.size()*sizeof(double);
    void *pointer = avx_data.data();
    align_offset = (double *) std::align(32, 32, pointer, size) - avx_data.data();

    //copy
    for (int i=0; i<Nm; i++) {
      avx_data[i + align_offset + avx_Nm*iAxi0] = Ak[i] * xi[i].x;
      avx_data[i + align_offset + avx_Nm*iAxi1] = Ak[i] * xi[i].y;
      avx_data[i + align_offset + avx_Nm*iAxi2] = Ak[i] * xi[i].z;

      // the cosine implementation computes cos(pi*x), so we'll divide out the pi here
      avx_data[i + align_offset + avx_Nm*ikkappa0] = k[i] / M_PI * kappa[i].x;
      avx_data[i + align_offset + avx_Nm*ikkappa1] = k[i] / M_PI * kappa[i].y;
      avx_data[i + align_offset + avx_Nm*ikkappa2] = k[i] / M_PI * kappa[i].z;

      // we also need to divide beta by pi, since that goes into the argument as well
      avx_data[i + align_offset + avx_Nm*ibeta] = beta[i] / M_PI;
    }
  }

  Vector3d TD13Field::getField(const Vector3d& pos) const {

#ifndef FAST_TD13
    Vector3d B(0.);
  
    for (int i=0; i<Nm; i++) {
      double z_ = pos.dot(kappa[i]);
      B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
    }

    return B;

#else // FAST_TD13

    // Initialize accumulators
    //
    // There is one accumulator per component of the result vector.
    // Note that each accumulator contains four numbers. At the end of
    // the loop, each of these number will contain the sum of every
    // fourth wavemodes, starting at a different offset. In the end, all
    // of the accumulator's numbers are added together (using
    // hsum_double_avx), resulting in the total sum.

    __m128d acc0 = _mm_setzero_pd();
    __m128d acc1 = _mm_setzero_pd();
    __m128d acc2 = _mm_setzero_pd();

    // broadcast position into SSE registers
    __m128d pos0 = _mm_set1_pd(pos.x);
    __m128d pos1 = _mm_set1_pd(pos.y);
    __m128d pos2 = _mm_set1_pd(pos.z);

    for (int i=0; i<avx_Nm; i+=2) {

      // load data from memory into AVX registers
      __m128d Axi0 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*iAxi0);
      __m128d Axi1 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*iAxi1);
      __m128d Axi2 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*iAxi2);

      __m128d kkappa0 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*ikkappa0);
      __m128d kkappa1 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*ikkappa1);
      __m128d kkappa2 = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*ikkappa2);

      __m128d beta = _mm_load_pd(avx_data.data() + i + align_offset + avx_Nm*ibeta);


      // Do the computation

      // this is the scalar product between k*kappa and pos
      __m128d z = _mm_add_pd(_mm_mul_pd(pos0, kkappa0),
			       _mm_add_pd(_mm_mul_pd(pos1, kkappa1),
					     _mm_mul_pd(pos2, kkappa2)
					     )
			       );

      // here, the phase is added on. this is the argument of the cosine.
      __m128d cos_arg = _mm_add_pd(z, beta);

      // ********
      // * Computing the cosine
      // *
      // * argument reduction
      // step 1: compute round(x), and store it in q
      __m128d q = _mm_round_pd(cos_arg, (_MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC));

      // we now compute s, which will be the input parameter to our polynomial approximation
      // of cos(pi/2*x) between 0 and 1
      // the andnot_pd is just a fast way of taking the absolute value
      __m128d s = _mm_sub_pd(cos_arg, q);
    
      // the following based on the int extraction process described here:
      // https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx/41223013
      // we assume -2^51 <= q < 2^51 for this, which is unproblematic, as double precision
      // has decayed far enough at that point that the cosine would be useless anyway.
    
      // we now want to check whether q is even or odd, because the cosine is negative for odd qs, so we'll have to flip the final result.
      // on an int, this is as simple as checking the 0th bit.
      // => manipulate the double in such a way that we can do this.
      // so, we add 2^52, such that the last digit of the mantissa is actually in the ones position.
      // since q may be negative, we'll also add 2^51 to make sure it's positive.
      // note that 2^51 is even and thus leaves evenness invariant, which is the only thing we care about here.

      q = _mm_add_pd(q, _mm_set1_pd(0x0018000000000000));

      // unfortunately, integer comparisons were only introduced in avx2, so we'll have to make do
      // with a floating point comparison to check whether the last bit is set.
      // however, masking out all but the last bit will result in a denormal float,
      // which may either result in performance problems or just be rounded down to zero,
      // neither of which is what we want here. To fix this, we'll mask in not only bit 0,
      // but also the exponent (and sign, but that doesn't matter) of q. Luckily, the exponent of q
      // is guaranteed to have the fixed value of 1075 (corresponding to 2^52) after our addition.

      __m128d invert = _mm_and_pd(q, _mm_castsi128_pd(_mm_set1_epi64x(0xfff0000000000001)));

      // if we did have a one in bit 0, our result will be equal to 2^52 + 1
      invert = _mm_cmpeq_pd(invert, _mm_castsi128_pd(_mm_set1_epi64x(0x4330000000000001)));

      // finally, we need to turn invert and right_invert into masks for the sign bit on each final double, ie
      invert = _mm_and_pd(invert, _mm_set1_pd(-0.0));

      // TODO: clamp floats between 0 and 1? This would ensure that we never see inf's, but maybe we want that,
      // so that things dont just fail silently...

      // * end of argument reduction
      // *******


      // ******
      // * evaluate the cosine using a polynomial approximation
      // * the coefficients for this were generated using sleefs gencoef.c
      // * I have no idea what I'm doing, so these coefficients are probably far from optimal.
      // * However, they should be sufficient for this case.
      s = _mm_mul_pd(s, s);

      __m128d u = _mm_set1_pd(+0.2211852080653743946e+0);

      u = _mm_add_pd(_mm_mul_pd(u, s), _mm_set1_pd(-0.1332560668688523853e+1 ));
      u = _mm_add_pd(_mm_mul_pd(u, s), _mm_set1_pd(+0.4058509506474178075e+1 ));
      u = _mm_add_pd(_mm_mul_pd(u, s), _mm_set1_pd(-0.4934797516664651162e+1));
      u = _mm_add_pd(_mm_mul_pd(u, s), _mm_set1_pd(1.));
    
      // then, flip the sign of each double for which invert is not zero. since invert
      // has only zero bits except for a possible one in bit 63, we can xor it onto
      // our result to selectively invert the 63st (sign) bit in each double where invert is set.
      u = _mm_xor_pd(u, invert);

      // * end computation of cosine
      // **********

      // Finally, Ak*xi is multiplied on. Since this is a vector, the
      // multiplication needs to be done for each of the three
      // components, so it happens separately.
      acc0 = _mm_add_pd(_mm_mul_pd(u, Axi0), acc0);
      acc1 = _mm_add_pd(_mm_mul_pd(u, Axi1), acc1);
      acc2 = _mm_add_pd(_mm_mul_pd(u, Axi2), acc2);
  }
  
  return Vector3d(hsum_double_sse(acc0),
                  hsum_double_sse(acc1),
                  hsum_double_sse(acc2)
                  );
#endif // FAST_TD13
}

} // namespace crpropa
