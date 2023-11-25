extern double erf();
double derf(x)
double *x;
{
  double derf;
  derf = erf(*x);
  return derf;
}
