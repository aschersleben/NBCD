#include <Rcpp.h>
using namespace Rcpp;

//////////////////////////////////////////
// Source for the memoziation approach: //
// http://stackoverflow.com/a/36973875  //
//////////////////////////////////////////

// Memoization structure to hold the hash map
struct mem_map{
  // Init static map
  static std::map<int, double> create_map()
  {
    std::map<int, double> m;
    m.clear();
    return m;
  }
  // Name map
  static std::map<int, double> memo;
};

// Instantiate class in global scope
std::map<int, double> mem_map::memo = mem_map::create_map();

// Reset map
//
// [[Rcpp::export]]
void clear_mem(){
  mem_map::memo.clear();
};

// Get map values.
//
// [[Rcpp::export]]
std::map<int, double> get_mem(){
  return mem_map::memo;
};



// Waiting Time Calculation
//
// See \code{\link{getWaitingTime}} for further information.
//
// @param p \code{[numeric]}\cr
//   Values have to add up to 1, will be checked. Otherwise returns an error.
//
// @param p1,p2,i1,imem
//   Parameters for recursion, do not change manually!
//
// [[Rcpp::export]]
double kMcpp(NumericVector p,
             NumericVector p1 = NumericVector::create(),
             NumericVector p2 = NumericVector::create(),
             IntegerVector i1 = IntegerVector::create(),
             int imem = 0) {

  int plen = p.size();
  int p1len = p1.size();
  int i1len = i1.size();

  // When starting, control if sum of "p" values is equal to 1
  if (p1len == 0) {
    double psabs = std::abs(std::accumulate(p.begin(), p.end(), 0.0) - 1.0);
    if (psabs > sqrt(std::numeric_limits< double >::epsilon())) {
      Rcpp::stop("Sum of p not equal to 1");
    }
    p1 = p;
    p1len = plen;
    mem_map::memo.clear();
  }

  // Look at hash map if already computed
  if(mem_map::memo.count(imem) > 0)
    return mem_map::memo[imem];

  if (i1len == 0) {
    i1 = IntegerVector(plen);
    for (int k = 0; k < p1len; k++) i1(k) = k + 1;
  }
  int p2len = p2.size();

  double out = 1;
  if (p1len == 1) {
    out = 1 / p1[0];
  } else {
    for (int i = 0; i < p1len; i++) {
      NumericVector pn = p1;
      pn.erase(i);
      NumericVector p2n (p2len + 1);
      for (int j = 0; j < p2len; j++) {
        p2n[j] = p2[j];
      }
      p2n[p2len] = p1[i];
      IntegerVector i1n = i1;
      i1n.erase(i);

      // i2: Name of cell in memory ("hash")
      int i2 = 0;
      int kl = p1len - 1;
      for (int k = 0; k < kl; k++) {
        i2 = i2 + i1n(k) * pow(2, i1n(k));
      }
      i2 = p1len * pow(10, 8) + i2;

      // recursion step:
      out = out + p1[i] * kMcpp(p, pn, p2n, i1n, i2);
    }

    if (p2len > 0) {
      double den = 1 - std::accumulate(p2.begin(), p2.end(), 0.0);
      out = out / den;
    }
  }
  mem_map::memo[imem] = out;
  return out;
};
