/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file ecm.h
 * \brief A library to implement Elliptic Curve Method
 *      to factorize a composite number.
 */

#include "mplib.h"

/**
 * \brief Number of different B values to try
 *      during the factorization
 *
 * If a B value does not produce a factor,
 * another one is tried until the number of
 * trials exceed B_THRESHOLD.
 */
#define B_THRESHOLD 1000
/**
 * \brief Number of different curves to try
 *      during the factorization
 *
 * If a curve does not produce a factor,
 * another one is tried until the number of
 * trials exceed CRV_THRESHOLD.
 */
#define CRV_THRESHOLD 20

/**
 * \brief Factorizes the given composite using Elliptic Curve Method
 * @param[out] d factor of n
 * @param[in] n number to be factorized
 * @param[in] nl number of digits of n in base \f$2^W\f$
 * @return an integer, positive when ECM succeeds, 0 otw
 */
int ecm(ui d, ui n, ui_t nl);
