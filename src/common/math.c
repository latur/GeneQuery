double pv(int a, int b, int c, int d, double * f)
{
    return exp(f[a + b] + f[c + d] + f[a + c] + f[b + d] - f[a + b + c + d] - f[a] - f[b] - f[c] - f[d]);
}

double right(int a, int b, int c, int d, double * f)
{
    long double psum = 0.0;
    while(1) {
        psum += pv(a, b, c, d, f);
        if (c == 0 || b == 0) break;
        a++; b--; c--; d++;
    }
    return psum;
}
