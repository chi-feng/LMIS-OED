#ifndef ExpectedInformationEstimator_h
#define ExpectedInformationEstimator_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <memory>
#include <Eigen/Core>

#include "Experiment.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::map;
using std::shared_ptr;
using std::string;
using std::tuple;
using std::vector;

class ExpectedInformationEstimator
{
  private:

    /// \brief Shared pointer to an instance of an experiment
    shared_ptr<Experiment> experiment;

    /// \brief Flag indicating whether the member data structures have been allocated
    bool allocated;

    /// \brief Flag indicating whether the member data structures relating to multiple importance sampling have been
    // allocated
    bool allocatedMIS;

  public:

    /// \brief Number of outer Monte Carlo Samples for estimating the expected information gain
    int N;

    /// \brief Number of Monte Carlo samples for estimating the \f$N\f$ marginal likelihoods \f$p(y^{(i)}|d)\f$
    int M1;

    /// \brief Number of Monte Carlo samples for estimating the \f$N\f$ conditional likelihoods
    // \f$p(y^{(i)}|\theta^{(i)},d)\f$
    int M2;

    /// \brief String indicating the type of distribution used for importance sampling, e.g. "MVT", "MVN", "EXACT" (when
    // available)
    string biasingDistributionType;

    /// \brief Degrees of freedom for multivariate-t biasing distributions
    double dof;

    /// \brief Size of nugget added to diagonal of estimated posterior covariances
    double nugget;

    /// \brief Maxmimum number of components in mixture distribution
    int maxComponents;

    /// \brief Minimum number of components in mixture distribution
    int minComponents;

    /// \brief Flag to toggle importance sampling (default = false)
    bool useIS;

    /// \brief Flag to toggle multiple importance sampling (default = false)
    bool useMIS;

    bool debugCondLikelihood;

    bool useMarginal;

    /// \brief Flag to toggle sorting of prior samples by prior likelihood (default = false)
    bool useSortHeuristic;

    bool useReverseLikelihood;

    bool useMinSampleDistance;

    /// \brief Flag to reuse prior samples between evaluations (default = false)
    bool reusePriorSamples;

    /// \brief Flag to use exact posterior if provided (default = false)
    bool useExactPosterior;

    /// \brief Flag to use bias correction
    bool useBiasCorrection;

    /// \brief Flag to toggle displaying a progress bar (default = false)
    bool showProgressBar;

    bool verbose;

  private:

    /// \brief Execution time in milliseconds
    double executionTime;

    /// \brief \f$N\times M_1\f$ matrix of log integrands for the marginal likelihood estimator where the \f$(i, j)\f$
    // element is the log of the \f$j\f$th integrand for estimating \f$p(y^{(i)}|d)\f$
    MatrixXd logf_ML;

    /// \brief \f$N\times M_1\f$ matrix of log importance weights for the marginal likelihood estimator where the \f$(i,
    // j)\f$ element is the log of the \f$j\f$th importance weight for estimating \f$p(y^{(i)}|d)\f$
    MatrixXd logw_ML;

    /// \brief \f$N\times M_2\f$ matrix of log integrands for the conditional likelihood estimator where the \f$(i, j)\f$
    // element is the log of the \f$j\f$th integrand for estimating \f$p(y^{(i)}|\theta^{(i)},d)\f$
    MatrixXd logf_CL;

    /// \brief \f$N\times M_2\f$ matrix of log importance weights for the conditional likelihood estimator where the
    // \f$(i, j)\f$ element is the log of the \f$j\f$th importance weight for estimating \f$p(y^{(i)}|\theta^{(i)},d)\f$
    MatrixXd logw_CL;

    /// \brief Length \f$N\f$ vector of all \f$N\f$ estimates of the log marginal likelihood \f$\log p(y^{(i)}|d)\f$
    VectorXd logML;

    /// \brief Length \f$N\f$ vector of all \f$N\f$ estimates of the log marginal likelihood \f$\log
    // p(y^{(i)}|\theta^{(i)},d)\f$
    VectorXd logCL;

    /// \brief Length \f$N\f$ vector of all \f$N\f$ estimates of the log marginal likelihood \f$\log p(y^{(i)}|d)\f$
    VectorXd logMLexact;

    /// \brief Length \f$N\f$ vector of all \f$N\f$ estimates of the log marginal likelihood \f$\log
    // p(y^{(i)}|\theta^{(i)},d)\f$
    VectorXd logCLexact;

    /// \brief \f$N\times(n_\theta+n_\eta)\f$ matrix where each row is a sample from the prior distribution
    MatrixXd X_prior;

    /// \brief \f$N\times n_y\f$ matrix where the \f$i\f$th row is the model evaluated at the \f$i\f$th sample from the
    // prior distribution
    MatrixXd G_prior;

    /// \brief \f$N\times n_y\f$ matrix where the \f$i\f$th row is the model evaluated at the \f$i\f$th sample from the
    // prior distribution with observation error
    MatrixXd Y_prior;

    /// \brief Collection of \f$N\f$ \f$M_1\times (n_\theta+n_\eta)\f$ matrices where the \f$i\f$th matrix contains the
    // sample set \f$\mathcal{X}_\text{ML}^{(i)}\sim q_\text{ML}^{(i)}\f$
    vector<MatrixXd> X_ML;

    /// \brief Collection of \f$N\f$ \f$M_1\times n_y\f$ matrices where the \f$i\f$th matrix contains the model evaluated
    // at the samples in \f$\mathcal{X}_\text{ML}^{(i)}\f$
    vector<MatrixXd> G_ML;

    /// \brief Collection of \f$N\f$ \f$M_2\times (n_\theta+n_\eta)\f$ matrices where the \f$i\f$th matrix contains the
    // sample set \f$\mathcal{X}_\text{ML}^{(i)}\sim q_\text{CL}^{(i)}\f$
    vector<MatrixXd> X_CL;

    /// \brief Collection of \f$N\f$ \f$M_2\times n_y\f$ matrices where the \f$i\f$th matrix contains the model evaluated
    // at the samples in \f$\mathcal{X}_\text{CL}^{(i)}\f$
    vector<MatrixXd> G_CL;

    /// \brief \f$N\times M_1\f$ matrix where element \f$(i,j)\f$ is the log prior density of the \f$j\f$th sample of
    // \f$\mathcal{X}_\text{ML}^{(i)}\f$
    MatrixXd ML_priorLogDensities;

    /// \brief \f$N\times M_1\f$ matrix where element \f$(i,j)\f$ is the log density \f$q_\text{ML}^{(i)}\f$ evaluated at
    // the \f$j\f$th sample of \f$\mathcal{X}_\text{ML}^{(i)}\f$
    MatrixXd ML_biasLogDensities;

    /// \brief Length \f$N\f$ vector of log prior densities at the \f$i\f$th sample of \f$\mathcal{X}_\text{prior}\f$
    VectorXd priorLogDensities;

    /// \brief \f$N\times N\f$ matrix of log likelihoods \f$p(y^{(i)}|\theta^{(j)},\eta^{(j)},d)\f$
    MatrixXd priorLogLikelihoods;

    vector<VectorXd> postMean;
    vector<MatrixXd> postCov;

    vector<VectorXd> condMean;
    vector<MatrixXd> condCov;

    /// \brief Collection of \f$N\f$ biasing distributions used to estimate the marginal likelihoods
    vector<shared_ptr<Distribution>> q_ML;

    /// \brief Collection of \f$N\f$ biasing distributions used to estimate the conditional likelihoods
    vector<shared_ptr<Distribution>> q_CL;

    /// \brief Map of tuples (a, b) that correspond to a vector with element i \f$\to\f$
    // \f$q_\text{ML}^{(a)}(\mathcal{X}_\text{ML}^{(b,i)}\f$
    map<tuple<int, int>, VectorXd> ML_cache;

    /// \brief \f$N\times N\f$ matrix of \f$q_\text{ML}^{(i)}(\mathcal{X}_\text{prior}^{(j)})\f$
    map<int, VectorXd> ML_priorCache;

    /// \brief Mixture indices for every iteration
    vector<vector<int>> mixtureIndices;

    /// \brief Effective sample size when estimating posterior moments for each iteration
    VectorXd ESS;

  public:

    ExpectedInformationEstimator(const shared_ptr<Experiment> experiment);
    void   PrintInfo(std::ostream& stream);
    void   WriteToFile(const std::string& prefix, const VectorXd& design);
    double Evaluate(const VectorXd& design);
    double GetExecutionTime();

  private:

    void   Allocate();
    void   AllocateMIS();
    void   SortPriorSamples();
    void   EvaluatePriorSamples(const VectorXd& design);
    double EvaluateBiasCorrectedEIG();
    double EvaluateNoIS(const VectorXd& design);
    double EvaluateIS(const VectorXd& design);
    double EvaluateMIS(const VectorXd& design);
};

inline std::ostream& operator<<(std::ostream& os, std::shared_ptr<ExpectedInformationEstimator> estimator)
{
  os << "-N " << estimator->N << std::endl;
  os << "-M1 " << estimator->M1 << std::endl;
  os << "-M2 " << estimator->M2 << std::endl;
  os << "-biasingDistributionType " << estimator->biasingDistributionType << std::endl;
  os << "-dof " << estimator->dof << std::endl;
  os << "-nugget " << estimator->nugget << std::endl;
  os << "-useExactPosterior " << estimator->useExactPosterior << std::endl;
  os << "-useIS " << estimator->useIS << std::endl;
  os << "-useMIS " << estimator->useMIS << std::endl;
  os << "-maxComponents " << estimator->maxComponents << std::endl;
  os << "-minComponents " << estimator->minComponents << std::endl;
  os << "-useSortHeuristic " << estimator->useSortHeuristic << std::endl;
  os << "-useReverseLikelihood " << estimator->useReverseLikelihood << std::endl;
  os << "-useMinSampleDistance " << estimator->useMinSampleDistance << std::endl;
  os << "-useBiasCorrection " << estimator->useBiasCorrection << std::endl;
  os << "-debugCondLikelihood " << estimator->debugCondLikelihood << std::endl;
  os << "-useMarginal " << estimator->useMarginal << std::endl;
  return os;
}

#endif // ifndef ExpectedInformationEstimator_h