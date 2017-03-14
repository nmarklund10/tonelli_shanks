#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <errno.h>
#include <climits>
#include <random>
#include <chrono>
#include <iostream>
#include "gmp.h"
#include <vector>
#include <algorithm>
#include <fstream>

int legendre_symbol(const mpz_t a, const mpz_t p) {
	mpz_t top;
	mpz_init_set(top, a);
	mpz_t bottom;
	mpz_init_set(bottom, p);
	mpz_t i;
	mpz_init(i);
	int result = 1;
	while (true) {
		if (mpz_cmp_si(top, 1) == 0) { // top == 1
			mpz_clear(i);
			mpz_clear(top);
			mpz_clear(bottom);
			return result;
		}
		else if (mpz_cmp_si(top, 2) == 0) { // top == 2
			mpz_set_si(i, 8);
			mpz_mod(i, bottom, i);
			if (mpz_cmp_si(i, 1) != 0 && mpz_cmp_si(i, 7) != 0) // i % 8 != 1 && i % 8 != 7
				result *= -1;
			mpz_clear(i);
			mpz_clear(top);
			mpz_clear(bottom);
			return result;
		}
		else if (mpz_cmp_si(top, 3) == 0 && mpz_cmp_si(bottom, 3) != 0) { // top == 3 && bottom != 3
			mpz_set_si(i, 12);
			mpz_mod(i, bottom, i);
			if (mpz_cmp_si(i, 1) == 0 || mpz_cmp_si(i, 11) == 0) {// i % 12 == 1 || i % 12 == 11
				mpz_clear(i);
				mpz_clear(top);
				mpz_clear(bottom);
				return result;
			}
			if (mpz_cmp_si(i, 5) == 0 || mpz_cmp_si(i, 7) == 0) {// i % 12 == 5 || i % 12 == 7
				mpz_clear(i);
				mpz_clear(top);
				mpz_clear(bottom);
				return -result;
			}
		}
		else if (mpz_cmp_si(top, 5) == 0 && mpz_cmp_si(bottom, 5) != 0) { // top == 5 && bottom != 5
			mpz_set_si(i, 5);
			mpz_mod(i, bottom, i);
			if (mpz_cmp_si(i, 1) != 0 && mpz_cmp_si(i, 4) != 0) // i % 5 != 1 && i % 5 != 4
				result *= -1;
			mpz_clear(i);
			mpz_clear(top);
			mpz_clear(bottom);
			return result;
		}
		mpz_mod(top, top, bottom); // top = top % bottom
		mpz_set_si(i, 1);
		mpz_and(i, top, i); // i = top & 1
		if (mpz_cmp_si(i, 1) == 0) { // top is odd
			mpz_set_si(i, 4);
			mpz_mod(i, top, i); // i = top % 4
			if (mpz_cmp_si(i, 3) == 0) {
				mpz_set_si(i, 4);
				mpz_mod(i, bottom, i); // i = bottom % 4
				if (mpz_cmp_si(i, 3) == 0)
					result *= -1;
			}
			mpz_set(i, top); // switch top and bottom
			mpz_set(top, bottom);
			mpz_set(bottom, i);
			mpz_mod(top, top, bottom);
		}
		mpz_set_si(i, 1);
		mpz_and(i, top, i);
		if (mpz_cmp_si(i, 0) == 0) { // top is even
			mpz_tdiv_q_2exp(top, top, 1); // top = top >> 1
			mpz_set_si(i, 2);
			result *= legendre_symbol(i, bottom);
		}
	}
}

void fast_exp_mod(mpz_t result, const mpz_t a, const mpz_t exp, const mpz_t p) {
	mpz_t exponent;
	mpz_init_set(exponent, exp);
	mpz_t base;
	mpz_init_set(base, a);
	mpz_set_si(result, 1);
	mpz_t comp;
	mpz_init(comp);
	mpz_t one;
	mpz_init_set_si(one, 1);
	while (mpz_cmp_si(exponent, 0) > 0) { // exponent > 0
		mpz_and(comp, exponent, one);
		if (mpz_cmp_si(comp, 1) == 0) { // exponent & 1 == 0
			mpz_mul(result, base, result);
			mpz_mod(result, result, p); // result = (base * result) % p
		}
		mpz_tdiv_q_2exp(exponent, exponent, 1); //exponent = exponent >> 1
		mpz_mul(base, base, base);
		mpz_mod(base, base, p); // base = (base * base) % p
	}
	mpz_clear(exponent);
	mpz_clear(base);
	mpz_clear(one);
	mpz_clear(comp);
}

bool miller_rabin(const mpz_t s, const mpz_t m, const mpz_t p) {
	if (mpz_cmp_si(p, 3) == 0) // p == 3
		return true;
	std::random_device rand_dev;
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	gmp_randseed_ui(state, rand_dev());
	mpz_t a;
	mpz_init(a);
	mpz_t b;
	mpz_init(b);
	
	mpz_t max;
	mpz_init(max);
	mpz_sub_ui(max, p, 2); // max = p - 2
	
	mpz_t p_min_1;
	mpz_init(p_min_1);
	mpz_sub_ui(p_min_1, p, 1);
	
	mpz_t j;
	mpz_init(j);
	mpz_t j_fin;
	mpz_init_set(j_fin, s);
	mpz_sub_ui(j_fin, j_fin, 1); // j_fin = s - 1
	
	for (int i = 0; i < 5; ++i) {	//less than 1% chance of false positive
		while (true) {
			mpz_urandomm(a, state, max); // a = random number in [0, p - 2]
			if (mpz_cmp_si(a, 2) >= 0) //try again if a < 2
				break;
		}
		fast_exp_mod(b, a, m, p);
		bool prob_prime = false;
		if (mpz_cmp_si(b, 1) != 0 && mpz_cmp(b, p_min_1) != 0) { // b != 1 && b != p - 1
			mpz_set_si(j, 0);
			while (mpz_cmp(j, j_fin) < 0) {
				mpz_mul(b, b, b);
				mpz_mod(b, b, p); // b = (b * b) % p
				if (mpz_cmp_si(b, 1) == 0) { // b == 1
					mpz_clear(a);
					mpz_clear(b);
					mpz_clear(max);
					mpz_clear(p_min_1);
					mpz_clear(j);
					mpz_clear(j_fin);
					return false;
				}
				if (mpz_cmp(b, p_min_1) == 0) { // b == p - 1
					prob_prime = true;
					break;
				}
				mpz_add_ui(j, j, 1); // j++
			}
			if (!prob_prime) {
				mpz_clear(a);
				mpz_clear(b);
				mpz_clear(max);
				mpz_clear(p_min_1);
				mpz_clear(j);
				mpz_clear(j_fin);
				return false;
			}
		}
	}
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(max);
	mpz_clear(p_min_1);
	mpz_clear(j);
	mpz_clear(j_fin);
	return true;
}

bool tonelli_shanks(mpz_t a, const mpz_t p, const mpz_t s, const mpz_t m, mpz_t i, mpz_t z, mpz_t e, mpz_t c, mpz_t x, mpz_t t, mpz_t temp, mpz_t b) {
	int num_of_iterations = 0;
	int max = mpz_get_ui(s) - 1;
	//Set e, c, x, t
	mpz_set(e, s);
	
	fast_exp_mod(c, z, m, p);
	
	mpz_add_ui(i, m, 1);
	mpz_tdiv_q_2exp(i, i, 1);
	fast_exp_mod(x, a, i, p);
	
	fast_exp_mod(t, a, m, p);
	
	// Iteration
	while(mpz_cmp_si(t, 1) != 0) {
		++num_of_iterations;
		// i = 1
		mpz_set_si(i, 1);
		// Find largest i such that t^(2^i) = 1 mod p
		while (true) {
			mpz_set_si(temp, 2);
			mpz_pow_ui(temp, temp, mpz_get_ui(i));
			fast_exp_mod(temp, t, temp, p);
			if (mpz_cmp_si(temp, 1) == 0)
				break;
			mpz_add_ui(i, i, 1);
		}
		//Compute b;
		mpz_sub(temp, e, i);
		mpz_sub_ui(temp, temp, 1);
		mpz_ui_pow_ui (temp, 2, mpz_get_ui(temp));
		fast_exp_mod(b, c, temp, p);
		//Update x, c, t, e
		mpz_mul(x, b, x);
		mpz_mod(x, x, p);
		mpz_mul(c, b, b);
		mpz_mod(c, c, p);
		mpz_mul(t, t, c);
		mpz_mod(t, t, p);
		mpz_set(e, i);
	}
	if (num_of_iterations == max)
		return true;
	return false;
}

int main(int argc, char** argv) {
	mpz_t p;
	mpz_init(p);
	mpz_t a;
	mpz_init(a);
	mpz_t m;
	mpz_init(m);
	mpz_t s;
	mpz_init(s);
	mpz_t i;
	mpz_init(i);
	mpz_t z;
	mpz_init(z);
	mpz_t e;
	mpz_init(e);
	mpz_t c;
	mpz_init(c);
	mpz_t x;
	mpz_init(x);
	mpz_t t;
	mpz_init(t);
	mpz_t b;
	mpz_init(b);
	mpz_t temp;
	mpz_init(temp);
	
	std::random_device rd;
    std::mt19937 gen(rd()); 
	
	std::ifstream infile("primes.txt");
	std::string prime;
	//number of times algorithm does max iterations
	std::vector<int> bins(10, 0);
	while (infile >> prime) {
		mpz_set_str(p, prime.c_str(), 10);
		
		mpz_sub_ui(m, p, 1);
		mpz_tdiv_q_2exp(m, m, 1); // m = (p - 1) >> 2
		mpz_set_ui(s, 1);
		while (true) {
			mpz_set_si(i, 1);
			mpz_and(i, m, i);
			if (mpz_cmp_si(i, 1) == 0) // m & 1 == 1 
				break;
			mpz_tdiv_q_2exp(m, m, 1); // m = m >> 1
			mpz_add_ui(s, s, 1); // s = s + 1
		}

		mpz_sub_ui(temp, p, 1);
		std::uniform_int_distribution<long long> dis(1, mpz_get_ui(temp));
		std::vector<long long> random_a;
		float iterations = 1000;
		for (int k = 0; k < iterations; ++k) {
			long long num = dis(gen);
			mpz_set_si(a, num);
			while(legendre_symbol(a, p) != 1) {
				num = dis(gen);
				mpz_set_si(a, num);
			}
			random_a.push_back(num);
		}
	
		mpz_set_ui(z, 7);
		int freq = 0;
		for (int j = 0; j < iterations; ++j) {
			mpz_set_si(a, random_a[j]);
			if (tonelli_shanks(a, p, s, m, i, z, e, c, x, t, temp, b))
				++freq;
		}
		
		if (freq < 100)
			++bins[0];
		else if (freq < 200)
			++bins[1];
		else if (freq < 300)
			++bins[2];
		else if (freq < 400)
			++bins[3];
		else if (freq < 500)
			++bins[4];
		else if (freq < 600)
			++bins[5];
		else if (freq < 700)
			++bins[6];
		else if (freq < 800)
			++bins[7];
		else if (freq < 900)
			++bins[8];
		else if (freq >= 900)
			++bins[9];		
	}
	
	for (size_t k = 0; k < 10; ++k) {
		std::cout << bins[k] << std::endl;
	}
	
	
	mpz_clear(b);
	mpz_clear(temp);
	mpz_clear(i);
	mpz_clear(z);
	mpz_clear(e);
	mpz_clear(c);
	mpz_clear(x);
	mpz_clear(t);
	mpz_clear(p);
	mpz_clear(a);
	mpz_clear(m);
	mpz_clear(s);
}

