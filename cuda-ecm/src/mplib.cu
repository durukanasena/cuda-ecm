/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file mplib.cu
 * \brief Implementation of mplib.h library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"

void big_rand(ui z, ui_t l) {
	int i;

	for(i = 0; i < l; i++) {
		z[i] = ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand());
	}
}

void big_mod_rand(ui z, ui_t l, ui n, ui_t nl, ui mu, ui_t mul) {
	ui_t z_[2 * nl];
	int i;

	for(i = 0; i < l; i++) {
		z_[i] = ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand()) * ((ui_t)rand());
	}
	for(i = l; i < 2 * nl; i++) {
		z_[i] = 0L;
	}
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_print(FILE *fp, ui a, ui_t al, char *s, char *R) {
    if(R != NULL) {
	    fprintf(fp, "%s := %s!(", s, R);
    } else {
        fprintf(fp, "%s := (", s);
    }
    fprintf(fp, "%u", a[0]);
    for(int i = 1; i < al; i++) {
        fprintf(fp, " + %u * (2^%d)^%d", a[i], W, i);
    }
    fprintf(fp, ");\n\n");
}

void big_is_equal(int *z, ui a, ui b, ui_t l) {
	int i;
	*z = 1;
	for(i = 0; i < l; i++) {
		if(a[i] != b[i]) {
			*z = 0;
		}
	}
}

void big_is_equal_ui(int *z, ui a, ui_t al, ui_t b) {
	int i;
	*z = 1;
	if(a[0] != b) {
		*z = 0;
	}
	for(i = 1; i < al; i++) {
		if(a[i] != 0) {
			*z = 0;
		}
	}
}

void big_add(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i;
	ui_t carry_bit = 0;

	for(i = 0; i < al; i++) {
		z[i] = a[i] + b[i] + carry_bit;
		if(z[i] < a[i]) {
			carry_bit = 1;
		} else if(z[i] > a[i]) {
			carry_bit = 0;
		}
	}
}

void big_mod_add(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu, ui_t mul) {
	int i;
	ui_t z_[2 * nl], carry_bit = 0;

	for(i = al + 1; i < 2 * nl; i++) {
		z_[i] = 0L;
	}
	for(i = 0; i < al; i++) {
		z_[i] = a[i] + b[i] + carry_bit;
		if(z_[i] < a[i]) {
			carry_bit = 1;
		} else if(z_[i] > a[i]) {
			carry_bit = 0;
		}
	}
	z_[al] = carry_bit;
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_sub(ui z, int *d, ui a, ui_t al, ui b, ui_t bl) {
	int i;
	ui_t borrow_bit = 0;

	for(i = 0; i < al; i++) {
		z[i] = a[i] - b[i] - borrow_bit;
		if(z[i] < a[i]) {
			borrow_bit = 0;
		} else if(z[i] > a[i]) {
			borrow_bit = 1;
		}
	}
	*d = borrow_bit;
}

void big_mod_sub(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl) {
	int i;
	ui_t z_[nl], borrow_bit = 0;

	for(i = 0; i < al; i++) {
		z_[i] = a[i] - b[i] - borrow_bit;
		if(z_[i] < a[i]) {
			borrow_bit = 0;
		} else if(z_[i] > a[i]) {
			borrow_bit = 1;
		}
	}
	if(borrow_bit) {
    	big_add(z, z_, nl, n, nl);
	} else {
		big_cpy(z, z_, 0, nl);
	}
}

void big_mul(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i, j;
	ui_t u, v;
	uni_t uv;

	for(i = 0; i <= al; i++) {
		z[i] = 0;
	}
 	for(i = 0; i < al; i++) {
	 	u = 0;
		for(j = 0; j < bl; j++) {
			uv = (uni_t)z[i + j] + (uni_t)a[i] * (uni_t)b[j] + (uni_t)u;
			u = uv >> W;
			v = uv & 0xFFFFFFFF;  // TODO: W != 32?
			z[i + j] = v;
		}
		z[i + bl] = u;
	}
}

void big_mod_mul(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu, ui_t mul) {
	int i, j;
	ui_t u, v, z_[2 * nl];
	uni_t uv;

	for(i = al + bl; i < 2 * nl; i++) {
		z_[i] = 0;
	}
	for(i = 0; i <= al; i++) {
		z_[i] = 0;
	}
 	for(i = 0; i < al; i++) {
	 	u = 0;
		for(j = 0; j < bl; j++) {
			uv = (uni_t)z_[i + j] + (uni_t)a[i] * (uni_t)b[j] + (uni_t)u;
			u = uv >> W;
			v = uv & 0xFFFFFFFF; // TODO: W != 32?
			z_[i + j] = v;
		}
		z_[i + bl] = u;
	}
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_get_mu(ui mu, ui n, ui_t nl) {
	mpz_t mp_n, mp_b2k, mp_mu;

    mpz_init(mp_n);
    mpz_init(mp_b2k);
    mpz_init(mp_mu);

	mpz_set_ui(mp_b2k, 0L);
	mpz_set_ui(mp_mu, 0L);

	mpz_import(mp_n, nl, -1, 4, 0, 0, n);
	mpz_add_ui(mp_b2k, mp_b2k, 1);
	mpz_mul_2exp(mp_b2k, mp_b2k, W * (2 * nl));
	mpz_fdiv_q(mp_mu, mp_b2k, mp_n);
	mpz_export(mu, NULL, -1, 4, 0, 0, mp_mu);
}

void big_get_A24(ui z, ui A, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
	ui_t c_2[nl], c_4[nl], A2[nl], ic_4[nl];
	int i, ret;

	c_2[0] = 2L;
	c_4[0] = 4L;
	for (i = 1; i < nl; i++) {
		c_2[i] = 0L;
		c_4[i] = 0L;
	}
	big_mod_add(A2, A, nl, c_2, nl, n, nl, mu, mul);
	ret = big_invert(ic_4, c_4, nl, n, nl);
	if(ret) { // Inverse exists
		big_mod_mul(z, A2, nl, ic_4, nl, n, nl, mu, mul);
		*flag = 1;
	} else { // Inverse does not exist
		big_gcd(z, nl, c_4, nl, n, nl);
		*flag = 0;
	}
}

uni_t barret_reduction_UL(uni_t p, uni_t b, uni_t k, uni_t z, uni_t m, uni_t L) { // Calculate z mod p where z < 2^W and p < 2^W
    uni_t bkpp = (k + 1) * L;
    uni_t bkmp = (k - 1) * L;
    uni_t bkp = 1 << bkpp;
    uni_t bkm = 1 << bkmp;

    uni_t q = ((z >> bkmp) * m) >> bkpp;
    uni_t r = (z & (bkp - 1)) - ((q * p) & (bkp - 1));
    if(r < 0) {
        r = r + bkp;
    }
    while(r >= p) {
        r -= p;
    }

    return r;
}

// ml = 2 * nl
void barret_reduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu, ui_t mul) { // Calculate m mod n
    ui_t k = nl, md[k + 1], mdmu[mul + k + 1], q[mul], mm[k + 1], qn[mul + nl], qnm[k + 1], r2[k + 1], r3[k + 1];
    int i, b;

    big_cpy(md, m, k - 1, k + 1); // md = m / b^(k - 1)
    big_mul(mdmu, md, k + 1, mu, mul); // mdmu = md * mu
    big_cpy(q, mdmu, k + 1, mul); // q = (m / b^(k - 1) * mu) / b^(k + 1)
    big_cpy(mm, m, 0, k + 1); // mm = m mod b^(k + 1)
    big_mul(qn, q, mul, n, nl); // qn = q * n
    big_cpy(qnm, qn, 0, k + 1); // qnm = qn mod b^(k + 1)
    big_sub(r3, &i, mm, k + 1, qnm, k + 1); // r3 = mm - qnm
	big_cpy(z, r3, 0, k);
    big_sub(r2, &b, r3, nl, n, nl); // while r >= n do: r <- r - n
	r2[nl] = r3[nl] - b;
    while(!(r2[nl] >> (W - 1))) {
        big_cpy(z, r2, 0, k);
		big_cpy(r3, r2, 0, k + 1);
        big_sub(r2, &b, r3, nl, n, nl);
		r2[nl] = r3[nl] - b;
    }
}

// Using GMP for now
void big_gcd(ui d, ui_t dl, ui a, ui_t al, ui b, ui_t bl) {
    mpz_t mp_a, mp_b, mp_d;
	int i;

	for(i = 0; i < dl; i++) {
		d[i] = 0L;
	}
	mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_d);

    mpz_set_ui(mp_a, 0L);
	mpz_set_ui(mp_b, 0L);
    mpz_set_ui(mp_d, 0L);

	mpz_import(mp_a, al, -1, 4, 0, 0, a);
    mpz_import(mp_b, bl, -1, 4, 0, 0, b);
	mpz_gcd(mp_d, mp_a, mp_b);
	mpz_export(d, NULL, -1, 4, 0, 0, mp_d);
}

int big_invert(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i, ret;
	mpz_t mp_z, mp_a, mp_b;

	mpz_init(mp_z);
	mpz_init(mp_a);
	mpz_init(mp_b);

	mpz_import(mp_a, al, -1, 4, 0, 0, a);
	mpz_import(mp_b, bl, -1, 4, 0, 0, b);
	mpz_set_ui(mp_z, 0L);

	for(i = 0; i < bl; i++) {
		z[i] = 0L;
	}
	ret = mpz_invert(mp_z, mp_a, mp_b);
	mpz_export(z, NULL, -1, 4, 0, 0, mp_z);       // iY2Z = Inv(Y^2Z)

	return ret;
}

void memoryAllocationGPU (ui *deviceArray, ui_t arraySize) {
	//Memory Allocation for GPU
	cudaMalloc(deviceArray, arraySize  * sizeof(ui_t));
	if (deviceArray == NULL)
		printf("deviceArray no space");
}

//__device__ void bigCpy(ui z, ui a, ui_t start, ui_t end, ui_t size) {
//	int i,j;
//	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
//	if(t_index < size){
//		int firstIndexZ = t_index * end, lastIndexZ = firstIndexZ + end;
//		int firstIndexA = t_index * start;
//		if(t_index == 0){
//			firstIndexA += start;
//		}
//		for(i = firstIndexZ, j = firstIndexA; i < lastIndexZ; i++, j++){
//			z[i] = a[j];
//		}
//	}
//}

__device__ void bigCpy(ui z, ui a, ui_t start, ui_t end, ui_t size) {
	int i,j;
	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
	if(t_index < size){
		int firstIndexZ = t_index * end, lastIndexZ = firstIndexZ + end;
		int firstIndexA = t_index * start;
		if(t_index == 0){
			firstIndexA += start;
		}
		bigCopy(z, firstIndexZ, lastIndexZ, a, firstIndexA);
	}
}

__device__ void bigMul(ui z, ui a, ui_t al, ui b, ui_t bl, ui_t size) {
	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
	int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
	ui_t u,low, high;
	if (t_index < size) {
		for(int i = firstIndexA; i <= lastIndexA; i++){
			z[i] = (ui_t)0;
		}
		for (int i = firstIndexA; i < lastIndexA; i++) {
			u = 0;
			for (int j = firstIndexA; j < lastIndexA; j++) {
				low = 0;
				high = 0;
				__mul_lo(low, a[i], b[j]);
				__mul_hi(high, a[i], b[j]);
				__add_cc(low, z[i+j], low);
				__addcy2(high);
				__add_cc(low, u, low);
				__addcy2(high);

				z[i+j] = low;
				u = high;
			}
			z[i + lastIndexA] = u;
		}
	}
}

__device__ void bigSub(ui z, ui controlBits, ui a, ui_t al, ui b , ui_t bl, ui_t size) {
	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
	int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
	int borrowbit = 0;

	if(t_index < size){
		for (int j = firstIndexA; j < lastIndexA; j++) {
			// c[j] = a[j] - b[j] - borrowBit
			__sub_cc(z[j], a[j], b[j]);
			__sub_cc(z[j], z[j], borrowbit);
			if(z[j]>a[j]){
				borrowbit = 1;
			}else if (z[j] < a[j]){
				borrowbit = 0;
			}
		}
		//__subcy(controlBits[t_index]);  //TODO : Add Macro
		controlBits[t_index] = borrowbit;
	}
}

__device__ float gpuGenerate(curandState* globalState) {
	int t_index = threadIdx.x + blockIdx.x * blockDim.x;
	curandState localState = globalState[t_index];
	float RANDOM = curand_uniform(&localState);
	globalState[t_index] = localState;
	return RANDOM;
}


__global__ void bigModRand(ui z, ui_t l, ui n, ui_t nl, ui mu, ui_t mul, ui_t size, curandState* globalState, ui_t seed){
	int t_index = threadIdx.x + blockIdx.x * blockDim.x;
	int i;

	if(t_index < size){
		ui z_ = new ui_t[2 * nl];
		int firstIndexZ_ = t_index * 2 * nl, lastIndexZ_ = firstIndexZ_ + (2 * nl);
		int halfIndexZ_ = ((lastIndexZ_ - firstIndexZ_ ) / 2 ) + firstIndexZ_;
		curand_init(seed, t_index, 0, &globalState[t_index]);
		for(i = firstIndexZ_; i < halfIndexZ_; i++) {
			z_[i] = (ui_t)(gpuGenerate(globalState) *  4294967295);
		}
		for(i = halfIndexZ_; i < lastIndexZ_; i++) {
			z_[i] = 0L;
		}
		barretReduction(z, z_, 2 * nl, n, nl, mu, mul,size);
	}
}

__global__ void bigModMul(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu, ui_t mul, ui_t size) {
	int i, j;
	ui z_ = new ui_t[2 * nl];
	ui_t u, low, high;
	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
	int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
	int firstIndexZ = t_index * 2 * nl, lastIndexZ = firstIndexZ + (2 * nl);
	if (t_index < size) {
		for(i = firstIndexZ; i < lastIndexZ; i++){
			z_[i] = 0;
		}
		for (i = firstIndexA; i < lastIndexA; i++) {
			u = 0;
			for (j = firstIndexA; j < lastIndexA; j++) {

				low = 0;
				high = 0;
				__mul_lo(low, a[i], b[j]);
				__mul_hi(high, a[i], b[j]);
				__add_cc(low, z_[i+j], low);
				__addcy2(high);
				__add_cc(low, u, low);
				__addcy2(high);

				z_[i + j] = low;
				u = high;
			}
			z_[i + lastIndexA] = u;
		}
		barretReduction(z, z_, 2 * nl, n, nl, mu, mul, size);
	}
}

__device__ void barretReduction(ui z, ui m, ui_t ml, ui n,ui_t nl, ui mu,ui_t mul, ui_t size){
    ui_t k = nl;
    ui md = new ui_t[k+1];
    ui mdmu = new ui_t[mul + k + 1];
    ui q = new ui_t[mul];
    ui mm = new ui_t[k + 1];
    ui qn = new ui_t[mul + nl];
    ui qnm = new ui_t[k + 1];
    ui r2 = new ui_t[k + 1];
    ui r3 = new ui_t[k + 1];
    //int i, j, b;
    ui controlBits = new ui_t[size];
    ui controlBitsi = new ui_t[size];
    int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

    bigCpy(md, m, k - 1, k + 1,size); // md = m / b^(k - 1)
    bigMul(mdmu, md, k + 1, mu, mul,size); // mdmu = md * mu
    bigCpy(q, mdmu, k + 1, mul,size); // q = (m / b^(k - 1) * mu) / b^(k + 1)
    bigCpy(mm, m, 0, k + 1,size); // mm = m mod b^(k + 1)
    bigMul(qn, q, mul, n, nl,size); // qn = q * n
    bigCpy(qnm, qn, 0, k + 1,size); // qnm = qn mod b^(k + 1)
    bigSub(r3, controlBitsi, mm, k + 1, qnm, k + 1,size); // r3 = mm - qnm
    bigCpy(z, r3, 0, k,size);
    bigSub(r2, controlBits, r3, nl, n, nl,size); // while r >= n do: r <- r - n
	r2[nl] = r3[nl] - controlBits[t_index];
    while(!(r2[nl] >> (W - 1))) {
    	bigCpy(z, r2, 0, k,size);
    	bigCpy(r3, r2, 0, k + 1,size);
        bigSub(r2, controlBits, r3, nl, n, nl,size);
		r2[nl] = r3[nl] - controlBits[t_index];
    }
}
