/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file test.h
 * \brief A library to test ecm
 */

/**
 * \brief Tests pro_curve_point function using GMP
 * 1. Generates a random curve and a projective point
 * using the function
 *
 * 2. Checks if the coefficients of the curve and the
 * point satisfies the curve equation using GMP
 */
void pro_curve_point_gmp_test();

/**
 * \brief Tests aff_curve_point function using GMP
 * 1. Generates a random curve and a affine point
 * using the function
 *
 * 2. Checks if the coefficients of the curve and the
 * point satisfies the curve equation using GMP
 */
void aff_curve_point_gmp_test();

/**
 * \brief Tests pro_add function using GMP
 * 1. Generates two random points
 *
 * 2. Computes the addition by implementing the
 * algebraic calculations using GMP
 *
 * 3. Compares the result of the function
 * with the GMP result
 */
void pro_add_gmp_test();

/**
 * \brief Tests pro_add function using Magma
 * 1.  Generates two random points
 *
 * 2. Computes the addition by implementing the
 * algebraic calculations using Magma
 *
 * 3. Compares the result of the function
 * with the Magma result
 */
void pro_add_magma_test();

/**
 * \brief Tests pro_dbl function using Magma
 * 1.  Generates a random point
 *
 * 2. Computes the double by implementing the
 * algebraic calculations using Magma
 *
 * 3. Compares the result of the function
 * with the Magma result
 */
void pro_dbl_magma_test();

// TODO: Implement maybe?
void pro_ladder_gmp_test();

// TODO: Implement
void pro_ladder_magma_test();

/**
 * \brief Tests ecm function
 * 1. Generates a random composite number
 *
 * 2. Calculates a factor of the number
 * using the function
 *
 * 3. Checks whether the factor found
 * actually divides the composite number
 */
void ecm_test();
