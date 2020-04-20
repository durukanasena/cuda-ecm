/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file montgomery.cu
 * \brief Implementation of montgomery.h library.
 */

#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"

#define THREADNUM 1
#define BLOCKNUM 1
#define SIZE 1

void pro_curve_point(ui d, MONTG_CURVE c, PRO_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A2[nl], X2[nl], AX2[nl], AX2Z[nl], X3[nl], Z2[nl], XZ2[nl], Y2[nl], Y2Z[nl], iY2Z[nl], RHS1[nl], RHS[nl];

    ui d_A2, d_X2, d_AX2, d_AX2Z, d_X3, d_Z2, d_XZ2, d_Y2, d_Y2Z, d_iY2Z, d_RHS1, d_RHS;
    ui d_A, d_B, d_X, d_Y, d_Z, d_n, d_mu;

    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Y = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);
    int i, s = 1, is_zero, is_one, is_n;

	memoryAllocationGPU(&d_A, nl);
	memoryAllocationGPU(&d_B, nl);
	memoryAllocationGPU(&d_X, nl);
	memoryAllocationGPU(&d_Y, nl);
	memoryAllocationGPU(&d_Z, nl);
	memoryAllocationGPU(&d_n, nl);
	memoryAllocationGPU(&d_mu, mul);

	memoryAllocationGPU(&d_A2, nl);
	memoryAllocationGPU(&d_X2, nl);
	memoryAllocationGPU(&d_AX2, nl);
	memoryAllocationGPU(&d_AX2Z, nl);
	memoryAllocationGPU(&d_X3, nl);
	memoryAllocationGPU(&d_Z2, nl);
	memoryAllocationGPU(&d_XZ2, nl);
	memoryAllocationGPU(&d_Y2, nl);
	memoryAllocationGPU(&d_Y2Z, nl);
	memoryAllocationGPU(&d_iY2Z, nl);
	memoryAllocationGPU(&d_RHS1, nl);
	memoryAllocationGPU(&d_RHS, nl);

	cudaMemset(d_A2, 0, nl * sizeof(ui_t));
	cudaMemset(d_X2, 0, nl * sizeof(ui_t));
	cudaMemset(d_AX2, 0, nl * sizeof(ui_t));
	cudaMemset(d_AX2Z, 0, nl * sizeof(ui_t));
	cudaMemset(d_X3, 0, nl * sizeof(ui_t));
	cudaMemset(d_Z2, 0, nl * sizeof(ui_t));
	cudaMemset(d_XZ2, 0, nl * sizeof(ui_t));
	cudaMemset(d_Y2, 0, nl * sizeof(ui_t));
	cudaMemset(d_Y2Z, 0, nl * sizeof(ui_t));
	cudaMemset(d_iY2Z, 0, nl * sizeof(ui_t));
	cudaMemset(d_RHS1, 0, nl * sizeof(ui_t));
	cudaMemset(d_RHS, 0, nl * sizeof(ui_t));

	srand(time(NULL));

    big_mod_rand(A, nl, n, nl, mu, mul);
    big_mod_rand(X, nl, n, nl, mu, mul);
    big_mod_rand(Y, nl, n, nl, mu, mul);
    big_mod_rand(Z, nl, n, nl, mu, mul);

//	cudaMemcpy(d_A, A, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_X, X, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_Y, Y, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
//	cudaMemcpy(d_Z, Z, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
//
//	curandState* globalStateA;
//	cudaMalloc (&globalStateA, nl * sizeof(curandState));
//	curandState* globalStateX;
//	cudaMalloc (&globalStateX, nl * sizeof(curandState));
//	curandState* globalStateY;
//	cudaMalloc (&globalStateY, nl * sizeof(curandState));
//	curandState* globalStateZ;
//	cudaMalloc (&globalStateZ, nl * sizeof(curandState));
//
//	bigModRand<<<BLOCKNUM,THREADNUM>>>(d_A, nl, d_n, nl, d_mu, mul, SIZE, globalStateA, rand());
//	bigModRand<<<BLOCKNUM,THREADNUM>>>(d_X, nl, d_n, nl, d_mu, mul, SIZE, globalStateX, rand());
//	bigModRand<<<BLOCKNUM,THREADNUM>>>(d_Y, nl, d_n, nl, d_mu, mul, SIZE, globalStateY, rand());
//	bigModRand<<<BLOCKNUM,THREADNUM>>>(d_Z, nl, d_n, nl, d_mu, mul, SIZE, globalStateZ, rand());
//
//	cudaMemcpy(A, d_A, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
//	cudaMemcpy(X, d_X, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
//	cudaMemcpy(Y, d_Y, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
//	cudaMemcpy(Z, d_Z, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
//	cudaMemcpy(n, d_n, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
//	cudaMemcpy(mu, d_mu, mul * sizeof(ui_t), cudaMemcpyDeviceToHost);


//    big_mod_mul(A2, A, nl, A, nl, n, nl, mu, mul);

    cudaMemcpy(d_A, A, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_A2, A2, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_n, n, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_mu, mu, mul * sizeof(ui_t), cudaMemcpyHostToDevice);

    bigModMul<<<BLOCKNUM,THREADNUM>>>(d_A2, d_A, nl, d_A, nl, d_n, nl, d_mu, mul, SIZE);

	cudaMemcpy(A2, d_A2, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(A, d_A, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(n, d_n, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(mu, d_mu, mul * sizeof(ui_t), cudaMemcpyDeviceToHost);

    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }                                                             // If A^2 mod n = 4, singular
    if(s != 0) {

    	cudaMemcpy(d_X, X, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
    	cudaMemcpy(d_X2, X2, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
    	cudaMemcpy(d_n, n, nl * sizeof(ui_t), cudaMemcpyHostToDevice);
    	cudaMemcpy(d_mu, mu, mul * sizeof(ui_t), cudaMemcpyHostToDevice);
    	bigModMul<<<BLOCKNUM,THREADNUM>>>(d_X2, d_X, nl, d_X, nl, d_n, nl, d_mu, mul, SIZE);
    	cudaMemcpy(X, d_X, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
    	cudaMemcpy(X2, d_X2, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
    	cudaMemcpy(n, d_n, nl * sizeof(ui_t), cudaMemcpyDeviceToHost);
    	cudaMemcpy(mu, d_mu, mul * sizeof(ui_t), cudaMemcpyDeviceToHost);

//        big_mod_mul(X2, X, nl, X, nl, n, nl, mu, mul);          // X2 = X^2
        big_mod_mul(AX2, A, nl, X2, nl, n, nl, mu, mul);        // AX = AX^2
        big_mod_mul(AX2Z, AX2, nl, Z, nl, n, nl, mu, mul);      // AX2Z = AX^2Z
        big_mod_mul(X3, X2, nl, X, nl, n, nl, mu, mul);         // X3 = X^3
        big_mod_mul(Z2, Z, nl, Z, nl, n, nl, mu, mul);          // Z2 = Z^2
        big_mod_mul(XZ2, X, nl, Z2, nl, n, nl, mu, mul);        // XZ2 = XZ^2
        big_mod_mul(Y2, Y, nl, Y, nl, n, nl, mu, mul);          // Y2 = Y^2
        big_mod_mul(Y2Z, Y2, nl, Z, nl, n, nl, mu, mul);        // Y2Z = Y^2Z

        big_mod_add(RHS1, X3, nl, AX2Z, nl, n, nl, mu, mul);    // RHS1 = X^3 + AX^2Z
        big_mod_add(RHS, RHS1, nl, XZ2, nl, n, nl, mu, mul);    // RHS = X^3 + AX^2Z + XZ^2
        big_gcd(d, nl, Y2Z, nl, n, nl);                         // d = GCD(Y^2Z, n)
        big_is_equal_ui(&is_one, d, nl, 1L);
        big_is_equal(&is_n, d, n, nl);
        if(is_one || is_n) {
            big_invert(iY2Z, Y2Z, nl, n, nl);                   // iY2Z = Inv(Y2Z)
            big_mod_mul(B, RHS, nl, iY2Z, nl, n, nl, mu, mul);  // B = RHS * Inv(Y^2Z)
            big_is_equal_ui(&is_zero, B, nl, 0L);               // If B mod n = 0, singular
            if(!is_zero) {
                c->A = A;
                c->B = B;
                c->n = n;
                p->X = X;
                p->Y = Y;
                p->Z = Z;
                *flag = 1;
            } else {
                *flag = -1;
            }
        } else {
            *flag = 0;
        }
    } else {
        *flag = -1;
    }
}

void aff_curve_point(ui d, MONTG_CURVE c, AFF_POINT p, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
    ui_t A2[nl], x2[nl], x3[nl], Ax2[nl], y2[nl], iy2[nl], rhs1[nl], rhs[nl];
    ui A = (ui)malloc(sizeof(ui_t) * nl);
    ui B = (ui)malloc(sizeof(ui_t) * nl);
    ui x = (ui)malloc(sizeof(ui_t) * nl);
    ui y = (ui)malloc(sizeof(ui_t) * nl);
    int i, s = 1, is_zero, is_one, is_n;

    big_mod_rand(A, nl, n, nl, mu, mul);
    big_mod_rand(x, nl, n, nl, mu, mul);
    big_mod_rand(y, nl, n, nl, mu, mul);

    big_mod_mul(A2, A, nl, A, nl, n, nl, mu, mul);
    if(A2[0] == 4L) {
        s = 0;
    }
    for(i = 1; i < nl; i++) {
        if(A2[i] != 0L) {
            s = 1;
        }
    }                                                           // Check singularity
    if(s != 0) {
        big_mod_mul(x2, x, nl, x, nl, n, nl, mu, mul);          // x2 = x^2
        big_mod_mul(Ax2, A, nl, x2, nl, n, nl, mu, mul);        // Ax2 = Ax^2
        big_mod_mul(x3, x2, nl, x, nl, n, nl, mu, mul);         // x3 = x^3
        big_mod_mul(y2, y, nl, y, nl, n, nl, mu, mul);          // y2 = y^2
        big_mod_add(rhs1, x3, nl, Ax2, nl, n, nl, mu, mul);              // rhs1 = x^3 + Ax^2
        big_mod_add(rhs, rhs1, nl, x, nl, n, nl, mu, mul);               // rhs = x^3 + Ax^2 + x
        big_gcd(d, nl, y2, nl, n, nl);                          // d = GCD(y^2, n)
        big_is_equal_ui(&is_one, d, nl, 1L);
        big_is_equal(&is_n, d, n, nl);
        if(is_one || is_n) {
            big_invert(iy2, y2, nl, n, nl);                     // iy2 = Inv(y^2)
            big_mod_mul(B, rhs, nl, iy2, nl, n, nl, mu, mul);   // B = rhs * Inv(y^2)
            big_is_equal_ui(&is_zero, B, nl, 0L);
            if(!is_zero) {
                c->A = A;
                c->B = B;
                c->n = n;
                p->x = x;
                p->y = y;
                *flag = 1;
            } else {
                *flag = -1;
            }
        } else {
            *flag = 0;
        }
    } else {
        *flag = -1;
    }
}

void pro_add(PRO_POINT p, PRO_POINT p1, PRO_POINT p2, PRO_POINT pd, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a[nl], b[nl],  c[nl], d[nl], da[nl], cb[nl], e[nl], f[nl], e2[nl], f2[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);

    big_mod_add(a, p2->X, nl, p2->Z, nl, n, nl, mu, mul);            // a = X2 + Z2
    big_mod_sub(b, p2->X, nl, p2->Z, nl, n, nl);            // b = X2 - Z2
    big_mod_add(c, p1->X, nl, p1->Z, nl, n, nl, mu, mul);            // c = X1 + Z1
    big_mod_sub(d, p1->X, nl, p1->Z, nl, n, nl);            // d = X1 - Z1
    big_mod_mul(da, d, nl, a, nl, n, nl, mu, mul);          // da = d * a
    big_mod_mul(cb, c, nl, b, nl, n, nl, mu, mul);          // cb = c * b
    big_mod_add(e, da, nl, cb, nl, n, nl, mu, mul);                  // e = da + cb
    big_mod_sub(f, da, nl, cb, nl, n, nl);                  // f = da - cb
    big_mod_mul(e2, e, nl, e, nl, n, nl, mu, mul);          // e2 = e^2
    big_mod_mul(f2, f, nl, f, nl, n, nl, mu, mul);          // f2 = f^2
    big_mod_mul(X, pd->Z, nl, e2, nl, n, nl, mu, mul);      // X = Zd * e2
    big_mod_mul(Z, pd->X, nl, f2, nl, n, nl, mu, mul);      // Z = Xd * f2

    p->X = X;
    p->Z = Z;
}

// Reuqires division
void aff_add(AFF_POINT p, AFF_POINT p1, AFF_POINT p2, ui A, ui B, ui n, ui_t nl, ui mu, ui_t mul) {
    // ui_t c[nl], c2[nl], c3[nl], d[nl], d2[nl], d3[nl];
    // ui x = (ui)malloc(sizeof(ui_t) * nl);
    // ui y = (ui)malloc(sizeof(ui_t) * nl);
    // int i, s;

    // big_mod_sub(c, p2->y, nl, p1->y, nl, n, nl);             // c = y2 - y1
    // big_mod_sub(d, p2->x, nl, p1->x, nl, n, nl);             // d = x2 - x1
    // big_mod_mul(c2, c, nl, c, nl, n, nl, mu, mul);           // c2 = c^2
    // big_mod_mul(d2, d, nl, d, nl, n, nl, mu, mul);           // d2 = d^2
    // big_mod_mul(c3, c2, nl, c, nl, n, nl, mu, mul);          // c3 = c^3
    // big_mod_mul(d3, d2, nl, d, nl, n, nl, mu, mul);          // d3 = d^3
}

void pro_dbl(PRO_POINT p, PRO_POINT p1, ui A24, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a [nl], a2[nl], b[nl], b2[nl], c[nl], d[nl], e[nl];
    ui X = (ui)malloc(sizeof(ui_t) * nl);
    ui Z = (ui)malloc(sizeof(ui_t) * nl);

    big_mod_add(a, p1->X, nl, p1->Z, nl, n, nl, mu, mul);            // a = X + Z
    big_mod_mul(a2, a, nl, a, nl, n, nl, mu, mul);          // a2 = a^2
    big_mod_sub(b, p1->X, nl, p1->Z, nl, n, nl);            // b = X - Z
    big_mod_mul(b2, b, nl, b, nl, n, nl, mu, mul);          // b2 = b^2
    big_mod_sub(c, a2, nl, b2, nl, n, nl);                  // c = a2 - b2
    big_mod_mul(X, a2, nl, b2, nl, n, nl, mu, mul);         // X = a2 * b2
    big_mod_mul(d, A24, nl, c, nl, n, nl, mu, mul);         // d = a24 * c
    big_mod_add(e, b2 ,nl, d, nl, n, nl, mu, mul);                   // e = b2 + d
    big_mod_mul(Z, c, nl, e, nl, n, nl, mu, mul);           // Z = c * e

    p->X = X;
    p->Z = Z;
}

// Reuqires division
void aff_dbl(ui x, ui z, ui x1, ui y1, ui A, ui B, ui n) {
}

void pro_ladder(PRO_POINT p, PRO_POINT p1, ui A24, ui k, ui_t kl, ui n, ui_t nl, ui mu, ui_t mul) {
    ui_t a, x;
    PRO_POINT R0 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT R1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT R0_ = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT R1_ = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    R0->X = (ui)malloc(sizeof(ui_t) * nl);
    R0->Y = (ui)malloc(sizeof(ui_t) * nl);
    R0->Z = (ui)malloc(sizeof(ui_t) * nl);
    R0_->X = (ui)malloc(sizeof(ui_t) * nl);
    R0_->Y = (ui)malloc(sizeof(ui_t) * nl);
    R0_->Z = (ui)malloc(sizeof(ui_t) * nl);
    R1_->X = (ui)malloc(sizeof(ui_t) * nl);
    R1_->Y = (ui)malloc(sizeof(ui_t) * nl);
    R1_->Z = (ui)malloc(sizeof(ui_t) * nl);
    p->X = (ui)malloc(sizeof(ui_t) * nl);
    p->Z = (ui)malloc(sizeof(ui_t) * nl);
    int i, j;

    big_cpy(R0->X, p1->X, 0, nl);
    big_cpy(R0->Z, p1->Z, 0, nl);
    pro_dbl(R1, p1, A24, n, nl, mu, mul);

    a = k[kl - 1];
    j = 0;
    while(a > 0) {
        a = a >> 1;
        j++;
    }                                                       // Find the index of the first 1
    for(i = j - 2; i >= 0; i--) {
        x = 1;
        x <<= i;
        big_cpy(R0_->X, R0->X, 0, nl);
        big_cpy(R0_->Z, R0->Z, 0, nl);
        big_cpy(R1_->X, R1->X, 0, nl);
        big_cpy(R1_->Z, R1->Z, 0, nl);
        if(!(k[kl - 1] & x)) {
            pro_dbl(R0, R0_, A24, n, nl, mu, mul);
            pro_add(R1, R0_, R1_, p1, n, nl, mu, mul);
        } else {
            pro_add(R0, R0_, R1_, p1, n, nl, mu, mul);
            pro_dbl(R1, R1_, A24, n, nl, mu, mul);
        }
    }
    for (i = kl - 2; i >= 0; i--) {
        for(j = W - 1; j >= 0; j--) {
            x = 1;
            x <<= i;
            big_cpy(R0_->X, R0->X, 0, nl);
            big_cpy(R0_->Z, R0->Z, 0, nl);
            big_cpy(R1_->X, R1->X, 0, nl);
            big_cpy(R1_->Z, R1->Z, 0, nl);
            if(!(k[kl - 1] & x)) {
                pro_dbl(R0, R0_, A24, n, nl, mu, mul);
                pro_add(R1, R0_, R1_, p1, n, nl, mu, mul);
            } else {
                pro_add(R0, R0_, R1_, p1, n, nl, mu, mul);
                pro_dbl(R1, R1_, A24, n, nl, mu, mul);
            }
        }
    }
    big_cpy(p->X, R0->X, 0, nl);
    big_cpy(p->Z, R0->Z, 0, nl);
}

void aff_ladder(ui x, ui y, ui x1, ui y1, ui k, ui n) {
}

int pro_is_on_curve(ui A, ui B, ui X, ui Y, ui Z, ui n) {
	return 0;
}

int aff_is_on_curve(ui A, ui B, ui x, ui y, ui n) {
	return 0;
}
