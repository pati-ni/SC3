#include <RcppArmadillo.h>
#include <set>
#include <vector>
#include <algorithm>
//#include <omp.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;

//' Updates pairwise global coverage matrix
//'
//' @param coverage full coverage matrix
//' @param subsample a vector that maps the subsample to the previous global sample
// [[Rcpp::export]]
void update_coverage_matrix(arma::mat& coverage, arma::rowvec& subsample)
{
  unsigned int subsample_size = subsample.n_elem;
  
  std::cerr << "I'm in! " << std::endl;
  //#pragma omp parallel for
  for(size_t i1 = 0; i1 < subsample_size; i1++)
  {
    // i2 is private, c++ hell yeah
    for (size_t i2 = i1; i2 < subsample_size; i2++)
    {
      
      unsigned int index1 = subsample(i1);
      unsigned int index2 = subsample(i2);
      std::cout << index1 <<"," << index2 << std::endl;
      coverage(index1,index2)++;
      coverage(index2, index1)++;
    }
  }
}


//' Merges the subsample consensus matrix of to the global one 
//' to the original matrix
//' 
//' @param subsample_consensus the consensus matrix produced from the SC3 sub system
//' @param global_consensus the global consensus_matrix, no shared access
//' @param subsample_map the subsampling matrix
// [[Rcpp::export]]
void subsample_merge(arma::mat& subsample_consensus, arma::mat& global_consensus, arma::rowvec& subsample_map)
{
  size_t subsample_size = subsample_consensus.n_rows;
  //#pragma omp parallel for
  for(size_t i1 = 0; i1 < subsample_size; i1++)
  {
    for (size_t i2 = i1; i2 < subsample_size; i2++)
    {
      size_t gl1 = subsample_map(i1);
      size_t gl2 = subsample_map(i2);
      // Random Access in large array, this might be slow
      // Can not do anything about it though :(
      global_consensus(gl1, gl2) += subsample_consensus(i1, i2);
      global_consensus(gl2, gl1) += subsample_consensus(i1, i2);// Whatever, it is a waste of memory anyways..
    }
  }
}



//' Consensus matrix computation
//' 
//' Computes consensus matrix given cluster labels
//' 
//' @param dat a matrix containing clustering solutions in columns
//' @param K number of clusters
// [[Rcpp::export]]
arma::mat consmx(const arma::mat& dat, int K) 
{
  using namespace std;
  typedef vector< set<int> > Membership;
  mat res = zeros(dat.n_rows, dat.n_rows);
 
  vector <Membership> clusters;
 
  // Build index
  for (size_t cm = 0; cm < dat.n_cols; cm++) 
  {
    Membership cluster(K , set<int>() );
    for (size_t i = 0; i < dat.n_rows; i++) 
    {
      cluster[dat(i,cm) - 1].insert(i);
    }
    
    // Push cluster lists to cluster membership container
    clusters.push_back(cluster);
  }
 
  // Build consensus matrix
  for (size_t i1 = 0; i1 < clusters.size(); i1++)
  {
    
    // Cluster Result
    Membership& cr1 = clusters[i1];
    for (size_t i2 = i1 + 1; i2 < clusters.size(); i2++)
    {
      // Comparing clustering result
      Membership& cr2 = clusters[i2];
      // Iterate through individual clusters in cr1
      for (Membership::const_iterator c1 = cr1.begin();  c1 != cr1.end(); ++c1)
      {
        
        for (Membership::const_iterator c2 = cr2.begin();  c2 != cr2.end(); ++c2)
        {
          vector <int> common_members;
          set_intersection(c1->begin(), c1->end(), c2->begin(), c2->end(), back_inserter(common_members));
          for (size_t m1 = 0; m1 < common_members.size(); m1++)
          {
            for (size_t m2 = m1 + 1; m2 < common_members.size(); m2++)
            {
              res(common_members[m1], common_members[m2])++;
              res(common_members[m2], common_members[m1])++;
            }
          }
        }
      }
    }
  }

  // Set diagonal back to one.. (Why? Nobody knows..)
  for ( size_t i = 0; i < dat.n_rows; i++)
  {
    // is this legit??
    res(i,i) = 1;
    // check if armadillo support assignment operator
    // it should not matter that much, if not we can continue..
  }
  
  // Results should be now consistent with the original implementation, just faster
  
  res /= dat.n_cols;
  return res;
}

//' Matrix left-multiplied by its transpose
//' 
//' Given matrix A, the procedure returns A'A.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat tmult(arma::mat x) {
    return(x.t()*x);
}

//' Converts the distance matrix to adjacency matrix
//' 
//' Given matrix A, the procedure returns a transformed matrix A'.
//' 
//' @param x Numeric matrix.
// [[Rcpp::export]]
arma::mat distance_to_adjacency_mat(arma::mat A) {
    return exp(A*A.max());
}


