/* How to compute a Fibonacci number fast.
   Step 1. Plan on using Matrix multiplication.  
        a. compute (special) matrix ** n  (raised to the power of n)
        b. extract the value from the resulting matrix.
        c. Print value, yay!

   Why this technique?
   The pedagogial example to compute Fibonacci numbers is thru recursion.
    Mainly to provide an introduction to 
     * recursion,
     * transforming between recursion and iteration
     * the costs of each form of implementation; ie the memory cost of recursion
       could be O(n), but should the same as iteration, O(1).

   The problem with those is that the computation cost is O(n).  That is lots of work,
   and grows fast.

   There is a faster way to do it O(log n), and we do it here.  It shows the
   power of knowing a bit more algebra and not misfocusing.

   A note on floating point precision.
   Floating point representation is used as it lets the program implementation
   be compact, clean.  They also let one compute very large Fibonacci numbers
   approximately, with no change in code.  It's also a nicety "modern" CPUs
   support floating point better -- larger number of bits, more registers, and
   more multiply/add units.  Beauty, eh?

   The 
      Bits in number    Where this loses accuracy in the lower digits
        64 bits          thru Fibonacci(77). (78 drops a lower bit)
        80 bits          thru Fibonacci(92). (94 drops a lower bit)
       106 bits 
       128 bits

    Different compilers provided different sizes of real numbers
      Bits in number    Where
        64 bits          Microsoft Visual C/C++
        80 bits          Most compilers on x86 (use /Qlong-double flag with the Intel C++ compiler)
       106 bits          PowerPC, SPARCv9
       128 bits          GNU C (4.3 and later) on x86, HP-UX, Solaris/SPARC

   If you want all the digits of love, you'll need to visit the satellite of
   multiprecision.
*/
#include <stdio.h>
#include <stdlib.h>

/** optional niceties */
#ifdef _MSC_VER
#include <float.h>
#pragma fenv_access (on)
#endif

/** optional niceties to get as much precisions as we can */
#if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 3))
#define real __float128
#else
#define real long double
#endif

static real fibonacci(long index);
int main(int argc, char* argv[])
{
   long index;
   real val;
   unsigned int control_word;

   // Check for a value to use as an index
   if (argc < 1)
      return 0;
   // Get the argument
   index = atol(argv[1]);

#ifdef _MSC_VER
   // Another nicety.
   // Tell the Microsoft Visual C/C++ environment to use as much precision
   // as it can.
   // Save the configuration of the FPU
   _controlfp_s(&control_word, 0, 0);
   // Set mantissa to 64 bits (80-bit floating point number internally)
   _controlfp_s(&control_word, _PC_64, MCW_PC);
#endif

   // Compute
   val = fibonacci(index);
   // print results
   printf("%Lf", val);

   return 0;
}

/* Step 2: Square the matrix faster.
    The matrix has specific form:
         | a b |
         | b c |

    We're going use that to simplify the calculations, and make 
    the process just that much faster.

    We are going to focus on squaring the matrix for now.  This is an important
    special case, as it is the most common form of multiplying the matrices.

    The resulting matrix looks like
         | x y |
         | y z |

     x = a**2 + b**2
     y = a*b + b*c 
     z = b**2 + c**2

     We'll compute b**2 only once
*/
static __inline void matrix_to2(real* mat)
{
    real a = mat[0];
    real b = mat[1];
    real b2= b*b;
    real c = mat[2];

    mat[0] = a*a + b2;
    mat[1] = a*b + b*c;
    mat[2] = b2  + c*c;
}


/* Step 3: Multiply the matrices
    Step 2 focused on multiplying a matrix against itself fast.  We're going
    multiply two matrices together now.  Both have the same symmetrical
    pattern.
         | a b |     | d e |
         | b c |  *  | e f |

    And the resulting matrix looks like
         | x y |
         | z w |

    This simplifies into 
       x = a*d + b*e
       y = a*e + b*f
       z = b*d + c*e
       w = b*e + c*f

    y and z are defined as being the same, but don't look the same.
    We're going to skip z.
*/
    
static __inline void matrix_multiply(real* out, real const* in1, real const* in2)
{
    real a= in1[0];
    real b= in1[1];
    real c= in1[2];

    real d= in2[0];
    real e= in2[1];
    real f= in2[2];

    real be = b*e;
    out[0] = a*d + be;
    out[1] = a*e + b*f;
    out[2] = be  + c*f;
}


/* Step 4: Raise a matrix to a power really fast
   part a.
   Computing the matrix raised to a power is easier than multiplying it n times.
   We're going to use the Egyptian powers method.  Which you already know most of,
   due to binary.  Let's look at 17.  Written in binary 17 is
          2**4 + 2**0
          (16  + 1)

   To raise matrix the power of 17 is pretty easy
         (m ** 16) * (m ** 1)

   16 is the same as 2**4, and 1 is the same as 2**0.  We aren't adding the
   numbers either, we're multipling.

   To raise a matrix to the power of 37 we would
         (m ** 32) * (m ** 4) * (m ** 1)

   Did you notice that the exponents are powers of two?  If we convert 37 to
   binary it is
         2**5 + 2**2 + 2**0
        ( 32  +   4  + 1   )

   Each of those "bits' in the binary form corresponds to the exponent of the
   power of 37.  What we're going to do is compute the powers of two of the matrix
   them combine them.  This greatly reduces the multiplies.

   part b.
   And to raise a matrix to a power of 2, it is just a matter of squaring the
   matrix of the next lower power of two.
       m ** 2 =  m * m
       m ** 4 = (m**2) * (m ** 2)
       and so on.

   part c. 
   The special matrix that is used to start off the computation of fibonacci
   numbers is
        | 1 1 |
        | 1 0 |

   part d.
   The resulting pattern of the matrix, which we'll use to get the result of
   the computation, is
        | F(n+1)  F(n)  |
        | F(n)    F(n-1)|
   We're going to combine these two features now
*/
static real fibonacci(long index)
{
   real mat[3] = {1, 1, 0};
   real ret[3] = {1, 1, 0};
   
   // Handle the case of 0 
   if (index == 0)
      return 0.0;
   if (index == 1)
      return 1.0;
   index --;

   // Go thru and multiply against the results powers of two.. if the corresponding
   // bit is set in the binary form of the number.
   // We do this from least-significant to most significant, as it lets us form the
   // the power of two as we go along, by squaring the base matrix.
   for(;;)
   {
      // If this power of 2 (bit) is set, multiply the matrix by it
      if (index & 1)
      {
         matrix_multiply(ret, ret, mat);
      }
      // See if we are done, and any further work is pointless
      index >>= 1;
      if (!index)
         break;
      // Raise the matrix to a power of two
      matrix_to2(mat);
   }

   return ret[1];
}

      
