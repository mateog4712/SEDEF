
#include <boost/math/tools/roots.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <sstream>
#include <unordered_map>

using namespace std;
extern "C" double tau_c(double edit_error, int kmer_size, int MAX_ERROR,int MAX_EDIT_ERROR)
{
	const double ERROR_RATIO = (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR;
   double gap_error = std::min(1.0, ERROR_RATIO * edit_error);
   double a = (1 - gap_error) / (1 + gap_error);
   double b = 1 / (2 * std::exp(kmer_size * edit_error) - 1);
   return a * b;
}

extern "C" double solve_inverse_jaccard_c(int j, int kmer_size, int MAX_ERROR,int MAX_EDIT_ERROR)
{
	if (j == 0) return 1;
	if (j == 1) return 0;
	return boost::math::tools::newton_raphson_iterate([j, kmer_size,MAX_ERROR,MAX_EDIT_ERROR](double d){
		const double ERROR_RATIO = (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR;
		double E = exp(d * kmer_size);
		return make_tuple(
			((1 - d * ERROR_RATIO) / (1 + d * ERROR_RATIO)) * (1.0 / (2 * E - 1)) - j,
			2 * (- kmer_size * E + ERROR_RATIO - 2 * ERROR_RATIO * E + E * kmer_size * pow(d * ERROR_RATIO, 2)) /
				pow((2 * E - 1) * (1 + d * ERROR_RATIO), 2)
		);
	}, 0.10, 0.0, 1.0, numeric_limits<double>::digits);
}


extern "C" int relaxed_jaccard_estimate_c(int s, int kmer_size, int MAX_EDIT_ERROR,int MAX_ERROR, double result)
{
	
	if (result != -1) return result;

	using namespace boost::math;
	const double CI = 0.75;
	const double Q2 = (1.0 - CI) / 2; // one side interval probability

	result = ceil(s * tau_c(MAX_EDIT_ERROR, kmer_size,MAX_ERROR,MAX_EDIT_ERROR));
	for (; result >= 0; result--) {
		double d = solve_inverse_jaccard_c(result / s, kmer_size,MAX_ERROR,MAX_EDIT_ERROR); // returns edit error
		double x = quantile(complement(binomial(s, tau_c(d, kmer_size,MAX_ERROR,MAX_EDIT_ERROR)), Q2)); // inverse binomial 
		double low_d = solve_inverse_jaccard_c(x / s, kmer_size,MAX_ERROR,MAX_EDIT_ERROR);
		if (100 * (1 - low_d) < MAX_EDIT_ERROR) {
			result++; 
			break;
		}
	}
	result = max(result, 0.0);
	return result;
}


extern "C" int fun1(int a){
    return a*4;
}



// #include <boost/icl/interval_map.hpp>
// #include <boost/icl/interval_set.hpp>
// #include <iostream>

// typedef boost::icl::discrete_interval<int> Interval;
// typedef boost::icl::interval_map<int, std::set<std::pair<Interval, Interval>>> Subtree;
// typedef boost::icl::interval_map<int, Subtree> Tree;

// int main(){

//     Tree tree;
//     auto a = Interval(1, 10);
//     auto b = Interval(5, 50);
//     tree += make_pair(a, Subtree({b, {make_pair(a, b)}})) ;
//      a = Interval(2, 12);
//      b = Interval(5, 16);
//     tree += make_pair(a, Subtree({b, {make_pair(a, b)}})) ;
//     auto pf = tree.find(6);
//     if (pf == tree.end()) std::cout << "OK";
//     else std::cout << (*pf).first << std::endl;
//     for(auto i = (*pf).second.begin();i!=(*pf).second.end();i++){
//         auto k = (*i).second;
//         for (auto it2: k)
//             std::cout << (it2).first << "\t" << (it2).second << (*i).first <<std::endl;
//     }
//     return 0;
// }
// auto ittree.find(3)
