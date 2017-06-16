#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include "dtwclustpp.h"

namespace dtwclust {

class SparseDistmatIndices {
public:
    SparseDistmatIndices(int num_rows) : _num_rows(num_rows) {}

    Rcpp::IntegerMatrix getNewIndices(const Rcpp::IntegerVector& i,
                                      const Rcpp::IntegerVector& j,
                                      bool symmetric)
    {
        std::vector<int> new_i, new_j;

        for (int ii = 0; ii < i.length(); ii++) {
            for (int jj = 0; jj < j.length(); jj++) {
                int this_i = i[ii];
                int this_j = j[jj];
                if (symmetric && this_j > this_i) {
                    int swap = this_i;
                    this_i = this_j;
                    this_j = swap;
                }
                int ij = this_i + (this_j - 1) * _num_rows;
                bool new_index = existing_indices.find(ij) == existing_indices.end();
                if (new_index) {
                    existing_indices.insert(ij);
                    new_i.push_back(this_i);
                    new_j.push_back(this_j);
                }
            }
        }

        Rcpp::IntegerMatrix new_ids(new_i.size(), 2);
        for (int k = 0; k < new_i.size(); k++) {
            new_ids(k, 0) = new_i[k];
            new_ids(k, 1) = new_j[k];
        }

        return new_ids;
    }

    Rcpp::IntegerVector getExistingIndices() {
        Rcpp::IntegerVector existing_ids(existing_indices.size());
        int i = 0;
        for (int id : existing_indices) {
            existing_ids[i] = id;
            i++;
        }
        return existing_ids;
    }

private:
    int _num_rows;
    std::unordered_set<int> existing_indices;
};

RCPP_MODULE(SparseDistmatIndices) {
    using namespace Rcpp;

    class_<SparseDistmatIndices>("SparseDistmatIndices")

    .constructor<int>()

    .method("getNewIndices", &SparseDistmatIndices::getNewIndices)
    .method("getExistingIndices", &SparseDistmatIndices::getExistingIndices)
    ;
}

} // namespace dtwclust
