# JAGWAS
Joint Analysis of multiple phenotype GWAS

<br />
Current version: 1.0.0 

## Quick Installation 

Option 1: Use the binary executable file for Linux
* Download the binary file from the main directory or run the following line of code:
```
git clone https://github.com/baihongguo/JAGWAS
cd JAGWAS
chmod a+x JAGWAS
```

Option 2: Build JAGWAS Library Dependencies  
   * C++11 compiler or later 
   * BLAS/LAPACK. For Intel processors, we recommend that GEM be compiled with an optimized math routine library such as the Intel oneAPI Math Kernal Library to replace BLAS/LAPACK for optimal performance.  
   * Armadillo Library. 

<br />

To install JAGWAS, run the following lines of code:
 ```
 git clone https://github.com/baihongguo/JAGWAS
 cd JAGWAS
 cd src/  
 make  
 ```
<br />
<br />
<br />


## Dependencies
C/C++ Compiler
 * A compiler with C++11 (or later) support is required.
 
LAPACK and BLAS
 * The LAPACK (Linear Algebra PACKage) and BLAS (Basic Linear Algebra Subprograms) libraries are used for matrix operations in GEM.

Intel processors:
 * We recommend linking GEM to the Intel oneAPI Math Kernal Library (oneMKL), instead of classical BLAS/LAPACK, for a greater performance boost. This can be done by replacing -llapack and -lblas in the makefile with -lmkl_gf_lp64 -lmkl_sequential -lmkl_core before compiling.
  * It is important to compile with -lmkl_sequential since GEM already does multi-threading across SNPs.

AMD processors:
 * For AMD processors, OpenBLAS (-lopenblas) may be a better alternative.

Armadillo Library
 * The Armadillo library is used for matrix operations.


## Usage

### Running GEM

1. [Command Line Options](#command-line-options)  
1. [Input Files](#input-files)
1. [Output File Format](#output-file-format)
1. [Examples](#examples)

<br /> 
<br />

### Command Line Options

Once JAGWAS is installed, the executable can be used to run the program.  
