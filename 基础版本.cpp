#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <chrono>
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

const int64 limit_inv_1 = 12286;
const int64 limit_inv_2 = 40951;
const int64 limit_inv_3 = 65521;
const int64 limit_inv_4 = 114661;

const int64 deg = 2048;
int64 res;
int64 limit = 1; //
int64 L;         // 二进制的位数
int64 RR[2 * deg];
int64 a[deg], b[deg];
int64 a_in[2 * deg], b_in[2 * deg];
int64 ret_1[deg], ret_2[deg], ret_3[deg], ret_4[deg];
int64 wn_rec[4097];

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

void NTT(int64* A, int64 type, int64 G, int64 mod)
{
    for (int i = 0; i < limit; ++i)
        if (i < RR[i])
            swap(A[i], A[RR[i]]);
    for (int mid = 1; mid < limit; mid <<= 1)
    { // 原根代替单位根
        int64 wn = modpow(G, (mod - 1) / (mid * 2), mod);
        int64 step = 2048 / mid;
        if (type == -1)
            wn = modpow(wn, mod - 2, mod);

        for (int len = mid << 1, pos = 0; pos < limit; pos += len)
        {
            int64 w = 1;
            for (int k = 0; k < mid; k++, w = (w* wn) % mod)
            {
                if (type == 1)
                    w = wn_rec[step * k];
                else
                {
                    w = wn_rec[step * (2 * mid - k)];
                    int dg = 0;
                }
                int64 x = A[pos + k];
                int64 y = w * A[pos + mid + k] % mod;
                A[pos + k] = (x + y);// % mod;
                A[pos + k + mid] = ((x - y));// +mod) % mod;// (x - y + mod) % mod;
            }
        }
    }
    if (type == -1)
    {
        int64 limit_inv = 0;
        if (mod == m1)
            limit_inv = limit_inv_1;
        if (mod == m2)
            limit_inv = limit_inv_2;
        if (mod == m3)
            limit_inv = limit_inv_3;
        if (mod == m4)
            limit_inv = limit_inv_4;
        for (int i = 0; i < limit; ++i)
            A[i] = (A[i] * limit_inv);//% mod;
    }
}
void poly_mul(int64* ret, int64* a, int64* b, int64 deg1, int64 G, int64 mod)
{
    memset(a_in, 0, sizeof(a_in));
    memset(b_in, 0, sizeof(b_in));
    memcpy(a_in, a, deg * sizeof(int64));
    memcpy(b_in, b, deg * sizeof(int64));
    for (limit = 1, L = 0; limit <= deg1; limit <<= 1)
        L++;
    for (int i = 0; i < limit; ++i)
    {
        RR[i] = (RR[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    //初始化根数组
    for (int i = 0; i < 4097; i++)
    {
        wn_rec[i] = modpow(G, (mod - 1) * i / 4096, mod);
    }

    NTT(a_in, 1, G, mod);
    NTT(b_in, 1, G, mod);
    for (int i = 0; i < limit; ++i)
        a_in[i] = a_in[i] * b_in[i]% mod;
    NTT(a_in, -1, G, mod);
    
    for (int i = 0; i < deg; i++)
    {
        ret[i] = (a_in[i] - a_in[deg + i] + mod * 10000000) % mod;
    }
}

int main()
{
    for (int i = 0; i < deg; ++i)
    {
        a[i] = 1;
        b[i] = 1;
    }
    auto start = high_resolution_clock::now();
    for (int i = 0; i < 1000; i++)
    {
        poly_mul(ret_1, a, b, 2 * deg - 2, G, m1);
        poly_mul(ret_2, a, b, 2 * deg - 2, G, m2);
        poly_mul(ret_3, a, b, 2 * deg - 2, G, m3);
        poly_mul(ret_4, a, b, 2 * deg - 2, G, m4);
    }
    auto end = high_resolution_clock::now();
    cout << "用时：" << duration_cast<milliseconds>(end - start).count() << "ms" << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m1 - ret_1[i];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m2 - ret_2[i];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m3 - ret_3[i];
    }
    cout << endl;
    for (int i = 0; i <= 10; ++i)
    {
        cout << setw(6) << m4 - ret_4[i];
    }
    cout << endl;
}
