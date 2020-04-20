/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file ecm.cu
 * \brief Implementation of ecm.h library.
 */

#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"
#include "ecm.h"

int ecm(ui d, ui n, ui_t nl) {
    ui A24 = (ui)malloc(sizeof(ui_t) * nl);
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t mu[nl + 1], B[nl];
    int i, j, coeff, is_one, is_n, flag;

    big_get_mu(mu, n, nl);
    for(i = 0; i <nl; i++) {
        d[i] = 0L;
    }
    for(i = 0; i < CRV_THRESHOLD; i++) {
        do {
            pro_curve_point(d, c, p1, n, nl, mu, nl + 1, &flag);
        } while(flag == -1);
        if(flag == 0) {
            return 1;
        } else {
            big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
            if(flag == 1) {
                coeff = 1;
                for(j = 1; j < 11; j++) {
                    coeff *= j;
                }
                B[0] = (ui_t)coeff;
                for(j = 1; j < nl; j++) {
                    B[j] = 0L;
                }
                // TODO: How to choose B?
                // TODO: How big is B?
                for(j = 0; j < B_THRESHOLD; j++) {
                    pro_ladder(p, p1, A24, B, nl, n, nl, mu, nl + 1);
                    if(p->Z == 0) {
                        big_gcd(d, nl, B, nl, n, nl);
                        big_is_equal_ui(&is_one, d, nl, 1L);
                        big_is_equal(&is_n, d, n, nl);
                        if(!is_one && !is_n) {
                            return 2;
                            // FIXME: Does not find in ladder!
                        }
                        coeff++;
                        // TODO: How to increment B?
                    }
                }
            } else {
                return 3;
                // FIXME: Does not find in (A+2)/4!
            }
        }
    }
    return 0;
}
