#include <Rcpp.h>
using namespace Rcpp;

// // [[Rcpp::depends(digest)]]
// // #include <pmurhashAPI.h>


/////////////////////////////////////////////////
// Source: http://stackoverflow.com/a/36973875 //
/////////////////////////////////////////////////


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



// bool isZero(double x) { return x == 0; }
//
// std::vector<double> remZ(std::vector<double> p) {
//   p.erase(std::remove_if(p.begin(), p.end(), isZero),
//           p.end());
//   return p;
// }



// kMcpp function
//
// @param p p
// @param p1 p1
// @param p2 p2
// @param i1 i1
// @param imem imem
// @export
// [[Rcpp::export]]
double kMcpp(NumericVector p,
             NumericVector p1 = NumericVector::create(),
             NumericVector p2 = NumericVector::create(),
             IntegerVector i1 = IntegerVector::create(),
             int imem = 0) {

  int plen = p.size();
  int p1len = p1.size();
  int i1len = i1.size();
  if (p1len == 0) {
    // p = remZ(p);
    // plen = p.size();
    double psabs = std::abs(std::accumulate(p.begin(), p.end(), 0.0) - 1.0);
    if (psabs > sqrt(std::numeric_limits< double >::epsilon())) {
      Rcpp::stop("Sum of p not equal to 1");
    }
    p1 = p;
    p1len = plen;
    mem_map::memo.clear();
  }

  // Look if already computed
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
      // i2 = p1len + i2 * p1len;
      i2 = p1len * pow(10, 8) + i2;
      //
      // int seed = 0;
      // std::string txt2 = std::to_string(i2);
      // const char* txt = txt2.c_str();
      // // const char txt[3];
      // // sprintf(txt, "%d", i2);
      // uint32_t nChar = strlen(txt);
      // unsigned int i2m = PMurHash32(seed, txt, nChar);

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


// /*** R
// k2cpp(rep(1/5, 5))
// k2.mem(rep(1/5, 5))
// kMcpp(rep(1/5, 5))
// */
