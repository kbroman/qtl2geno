// intercross QTLCross class (for HMM)
// this is a hack for Tony Long, to deal with the case of super-low coverange sequencing
// at any given polymorphism, you just observe a single read, "A" or "B"
// emission probs modified to by Pr(A|het) = Pr(B|het) = 1/2

#ifndef CROSS_F2SEQ_H
#define CROSS_F2SEQ_H

#include <Rcpp.h>
#include "cross.h"

class F2SEQ : public QTLCross
{
 public:
    F2SEQ(){
        crosstype = "f2seq";
        phase_known_crosstype = "f2seqpk";
    };
    ~F2SEQ(){};

    const bool check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const double init(const int true_gen,
                      const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);
    const double emit(const int obs_gen, const int true_gen, const double error_prob,
                      const Rcpp::IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const Rcpp::IntegerVector& cross_info);
    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const Rcpp::IntegerVector possible_gen(const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const int ngen(const bool is_x_chr);

    const Rcpp::NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const Rcpp::LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr);

    const Rcpp::NumericMatrix get_x_covar(const Rcpp::LogicalVector& is_female, const Rcpp::IntegerMatrix& cross_info);

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles, const bool is_x_chr);
};

#endif // CROSS_F2SEQ_H
