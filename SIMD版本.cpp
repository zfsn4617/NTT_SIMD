#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <immintrin.h>
using namespace std;
using namespace std::chrono;
typedef int64_t int64;
typedef unsigned int uint;

const int64 m1 = 12289;
const int64 m2 = 40961;
const int64 m3 = 65537;
const int64 m4 = 114689;
const int64 M1 = (int64)m2 * m3 * m4;
const int64 M2 = (int64)m1 * m3 * m4;
const int64 M3 = (int64)m1 * m2 * m4;
const int64 M4 = (int64)m1 * m2 * m3;
const int64 m1_ = 11711;
const int64 m2_ = 37278;
const int64 m3_ = 20721;
const int64 m4_ = 94134;
const int64 mod = m1 * m2 * m3 * m4;
const int64 G = 23;
const int64 Gi_1 = 6946;
const int64 Gi_2 = 21371;
const int64 Gi_3 = 45591;
const int64 Gi_4 = 9973;

int64 wn_rec[4097];
__m256i wn_rec_avx[4097];
__m256i m_avx;// = { m1,m2,m3,m4 };
__m256i limit_inv;// = { 12286,40951,65521,114661 };

const int64 deg = 2048;
int64 res;
int64 limit = 1; //
int64 L;         // 二进制的位数
int64 RR[2 * deg];
int64 a[deg], b[deg];
int64 a_in[2 * deg], b_in[2 * deg];
int64 ret_1[deg], ret_2[deg], ret_3[deg], ret_4[deg];

__m256i a_avx[deg];
__m256i b_avx[deg];
__m256i a_in_avx[2 * deg];
__m256i b_in_avx[2 * deg];
__m256i ret_avx[deg];

int64 modpow(int64 a, int64 b, int64 mod)
{
    int64 res = 1;
    while (b)
    {
        if (b & 1)
            res = res * a % mod;
        a = a * a % mod;
        b >>= 1;
    }
    return res % mod;
}

int64 inv(int64 x, int64 mod)
{
    return modpow(x, mod - 2, mod);
}

void NTT(__m256i* A, int64 type, int64 G, __m256i mod)
{
    for (int i = 0; i < limit; ++i)
        if (i < RR[i])
            swap(A[i], A[RR[i]]);
    for (int mid = 1; mid < limit; mid <<= 1)
    { // 原根代替单位根
        int64 step = 2048 / mid;
        for (int len = mid << 1, pos = 0; pos < limit; pos += len)
        {
            for (int k = 0; k < mid; k++)
            {
                __m256i w;
                if (type == 1)
                    w = wn_rec_avx[step * k];
                else
                    w = wn_rec_avx[step * (2 * mid - k)];
                __m256i x = A[pos + k];
                __m256i y;// = w * A[pos + mid + k] % mod;
                y = _mm256_mul_epi32(w, A[pos + mid + k]);
                y.m256i_i64[0] %= mod.m256i_i64[0];
                y.m256i_i64[1] %= mod.m256i_i64[1];
                y.m256i_i64[2] %= mod.m256i_i64[2];
                y.m256i_i64[3] %= mod.m256i_i64[3];
                A[pos + k] = _mm256_add_epi64(x, y);// (x + y);
                A[pos + k + mid] = _mm256_sub_epi64(x, y);// ((x - y));
            }
        }
    }
    if (type == -1)
    {
        //int64 limit_inv = inv(limit, mod);
        for (int i = 0; i < limit; ++i)
            A[i] = _mm256_mul_epi32(A[i], limit_inv);//(A[i] * limit_inv);//% mod;
    }
}
void poly_mul(__m256i* ret, __m256i* a, __m256i* b, int64 deg1, int64 G, __m256i mod)
{
    memset(a_in_avx, 0, sizeof(a_in_avx));
    memset(b_in_avx, 0, sizeof(b_in_avx));
    memcpy(a_in_avx, a, deg * sizeof(__m256i));
    memcpy(b_in_avx, b, deg * sizeof(__m256i));
    for (limit = 1, L = 0; limit <= deg1; limit <<= 1)
        L++;
    for (int i = 0; i < limit; ++i)
    {
        RR[i] = (RR[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }

    NTT(a_in_avx, 1, G, m_avx);
    NTT(b_in_avx, 1, G, m_avx);
    for (int i = 0; i < limit; ++i)
    {
        a_in_avx[i] = _mm256_mul_epi32(a_in_avx[i], b_in_avx[i]);
        a_in_avx[i].m256i_i64[0] %= mod.m256i_i64[0];
        a_in_avx[i].m256i_i64[1] %= mod.m256i_i64[1];
        a_in_avx[i].m256i_i64[2] %= mod.m256i_i64[2];
        a_in_avx[i].m256i_i64[3] %= mod.m256i_i64[3];
            ;// a_in[i] * b_in[i] % mod;
    }
    NTT(a_in_avx, -1, G, m_avx);

    __m256i mod_24 = _mm256_slli_epi64(mod, 24);

    for (int i = 0; i < deg; i++)
    {
        __m256i temp = _mm256_sub_epi64(a_in_avx[i], a_in_avx[deg + i]);
        __m256i temp2 = _mm256_add_epi64(temp,mod_24);
        ret[i].m256i_i64[0] = temp2.m256i_i64[0] % mod.m256i_i64[0];
        ret[i].m256i_i64[1] = temp2.m256i_i64[1] % mod.m256i_i64[1];
        ret[i].m256i_i64[2] = temp2.m256i_i64[2] % mod.m256i_i64[2];
        ret[i].m256i_i64[3] = temp2.m256i_i64[3] % mod.m256i_i64[3];
    }
}

int main()
{
    //初始化全局变量
    m_avx.m256i_i64[0] = m1;
    m_avx.m256i_i64[1] = m2;
    m_avx.m256i_i64[2] = m3;
    m_avx.m256i_i64[3] = m4;
    limit_inv.m256i_i64[0] = 12286;
    limit_inv.m256i_i64[1] = 40951;
    limit_inv.m256i_i64[2] = 65521;
    limit_inv.m256i_i64[3] = 114661;
    //初始化根数组
    for (int i = 0; i < 4097; i++)
    {
        wn_rec_avx[i].m256i_i64[0] = modpow(G, (m1 - 1) * i / 4096, m1);
        wn_rec_avx[i].m256i_i64[1] = modpow(G, (m2 - 1) * i / 4096, m2);
        wn_rec_avx[i].m256i_i64[2] = modpow(G, (m3 - 1) * i / 4096, m3);
        wn_rec_avx[i].m256i_i64[3] = modpow(G, (m4 - 1) * i / 4096, m4);
    }

    for (int i = 0; i < deg; ++i)
    {
        for (int j = 0; j < 4; j++)
        {
            a_avx[i].m256i_i64[j] = 1;
            b_avx[i].m256i_i64[j] = 1;
        }
    }
    auto start = high_resolution_clock::now();
    for (int i = 0; i < 1000; i++)
    {
        poly_mul(ret_avx, a_avx, b_avx, 2 * deg - 2, G, m_avx);
    }
    auto end = high_resolution_clock::now();
    cout << "用时：" << duration_cast<milliseconds>(end - start).count() << "ms" << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m1 - ret_avx[i].m256i_i64[0];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m2 - ret_avx[i].m256i_i64[1];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m3 - ret_avx[i].m256i_i64[2];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m4 - ret_avx[i].m256i_i64[3];
    }
    cout << endl;
}
