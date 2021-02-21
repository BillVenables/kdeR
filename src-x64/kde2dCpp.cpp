#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List kde2dCpp(NumericVector x,  NumericVector y,
              NumericVector xo, NumericVector yo,
              NumericVector bwv, double area)
{
  int nx = xo.size(), ny = yo.size(), n = x.size();
  NumericMatrix z(nx, ny);
  double* xd = x.begin();
  double* yd = y.begin();
  double* bw = bwv.begin();

  for(int k = 0; k < n; k++ )
  {
    NumericVector dx = dnorm(xo, xd[k], bw[0]),
                  dy = dnorm(yo, yd[k], bw[1]);
    for(int j = 0; j < ny; j++)
    {
      z(_, j) = z(_, j) + dx*dy(j);
    }
  }
  return List::create(_["x"] = xo,
                      _["y"] = yo,
                      _["z"] = z/(sum(z) * area));
}
