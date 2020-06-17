#include <map>
#include <algorithm>

#include <cstdio>

#include "Distribution.h"
#include "MultivariateNormal.h"
#include "MultivariateT.h"

#include "Utilities.h"
#include "Timer.h"

#include "ExpectedInformationEstimator.h"

using namespace std;
using namespace Eigen;
using namespace Utilities;

ExpectedInformationEstimator::ExpectedInformationEstimator(const shared_ptr<Experiment> experiment) : experiment(experiment)
{
  useIS                   = false;
  useMIS                  = false;
  useSortHeuristic        = false;
  useReverseLikelihood    = true;
  useMinSampleDistance    = false;
  useBiasCorrection       = false;
  useMarginal             = false;
  verbose                 = false;
  biasingDistributionType = "MVT";
  dof                     = 3;
  nugget                  = 1e-4;
  N                       = 100;
  M1                      = 100;
  M2                      = 100;
  allocated               = false;
  allocatedMIS            = false;
  reusePriorSamples       = false;
  useExactPosterior       = false;
  maxComponents           = pow(2, experiment->inputDim);
  minComponents           = 1;
  executionTime           = 0;
  debugCondLikelihood     = false;
}

double ExpectedInformationEstimator::GetExecutionTime()
{
  return executionTime;
}

void ExpectedInformationEstimator::WriteToFile(const std::string& prefix, const VectorXd& design)
{
  cout << "ExpectedInformationEstimator::WriteToFile" << endl;
  // compute exact marginal and conditional likelihood (if possible)
  if (experiment->hasPosterior) {
    for (int i = 0; i < N; ++i) {
      VectorXd postMean = experiment->PosteriorMean(Y_prior.row(i), design);
      MatrixXd postCov  = experiment->PosteriorCovariance(Y_prior.row(i), design);
      VectorXd condMean;
      MatrixXd condCov;
      tie(condMean, condCov) = MultivariateNormal::GetConditional(postMean, postCov, X_prior.row(i).head(experiment->parametersOfInterest));
      auto     q_ML = make_shared<MultivariateNormal>(postMean, postCov);
      auto     q_CL = make_shared<MultivariateNormal>(condMean, condCov);
      int n         = 5; // zero-variance estimator but for good measure
      VectorXd x;
      VectorXd lml = VectorXd::Zero(n);
      VectorXd lcl = VectorXd::Zero(n);
      for (int j = 0; j < n; j++) {
        x                                      = q_ML->GetSample();
        lml(j)                                 = experiment->LogLikelihoodCached(Y_prior.row(i), x, experiment->Evaluate(x, design), design) + experiment->PriorLogDensity(x) - q_ML->LogDensity(x);
        x                                      = X_prior.row(i);
        x.tail(experiment->nuisanceParameters) = q_CL->GetSample();
        lcl(j)                                 = experiment->LogLikelihoodCached(Y_prior.row(i), x, experiment->Evaluate(x, design), design) + experiment->CondPriorLogDensity(x) - q_CL->LogDensity(x.tail(experiment->nuisanceParameters));
      }
      logMLexact(i) = LogSumExp(lml) - log(n);
      logCLexact(i) = LogSumExp(lcl) - log(n);
    }
    WriteEigenAsciiFile(prefix + "logMLexact", logMLexact);
    WriteEigenAsciiFile(prefix + "logCLexact", logCLexact);
  }

  // WriteEigenBinaryFile(prefix + "X_prior", X_prior);
  // WriteEigenBinaryFile(prefix + "logf_ML", logf_ML);
  // WriteEigenBinaryFile(prefix + "logw_ML", logw_ML);
  // WriteEigenBinaryFile(prefix + "logf_CL", logf_CL);
  // WriteEigenBinaryFile(prefix + "logw_CL", logw_CL);

  WriteEigenAsciiFile(prefix + "logML", logML);
  WriteEigenAsciiFile(prefix + "logCL", logCL);
  WriteEigenAsciiFile(prefix + "priorLogDensities", priorLogDensities);
  // WriteEigenBinaryFile(prefix + "priorLogLikelihoods", priorLogLikelihoods);
  WriteEigenAsciiFile(prefix + "ESS", ESS);

  // compute customized ess

  VectorXd cESS_ML = VectorXd::Zero(N);
  VectorXd cESS_CL = VectorXd::Zero(N);
  for (int i = 0; i < N; ++i) {
    cESS_ML(i) = CustomizedEffectiveSampleSize(logf_ML.row(i).array().exp(), logw_ML.row(i).array().exp());
    cESS_CL(i) = CustomizedEffectiveSampleSize(logf_CL.row(i).array().exp(), logw_CL.row(i).array().exp());
  }
  WriteEigenAsciiFile(prefix + "cESS_ML", cESS_ML);
  WriteEigenAsciiFile(prefix + "cESS_CL", cESS_CL);

  /*
   *  if (useMIS || useIS) {
   *  IOFormat full(FullPrecision, DontAlignCols, " ", "\n", "", "", "", "");
   *  ofstream file;
   *  file.open(prefix + "postMean");
   *  for (int i = 0; i < N; ++i) {
   *   file << postMean[i].transpose().format(full) << endl;
   *  }
   *  file.close();
   *  file.open(prefix + "postCov");
   *  for (int i = 0; i < N; ++i) {
   *   file << postCov[i].format(full) << endl;
   *  }
   *  file.close();
   *  if (experiment->hasPosterior) {
   *   file.open(prefix + "postMeanExact");
   *   for (int i = 0; i < N; ++i) {
   *     file << experiment->PosteriorMean(Y_prior.row(i), design).transpose().format(full) << endl;
   *   }
   *   file.close();
   *   file.open(prefix + "postCovExact");
   *   for (int i = 0; i < N; ++i) {
   *     file << experiment->PosteriorCovariance(Y_prior.row(i), design).format(full) << endl;
   *   }
   *  }
   *  }
   */

  if (useMIS) {
    ofstream out;
    out.open(prefix + "mixtureIndices");
    for (int i = 0; i < N; ++i) {
      out << i;
      for (size_t m = 1; m < mixtureIndices[i].size(); ++m) {
        out << " " << mixtureIndices[i][m];
      }
      out << endl;
    }
    out.close();
  }
  if (verbose) {
    ofstream file;
    file.open(prefix + "X_ML", ios::out | ios::binary);
    int rows = N * M2;
    int cols = experiment->inputDim;
    file.write(reinterpret_cast<char *>(&rows), sizeof(rows));
    file.write(reinterpret_cast<char *>(&cols), sizeof(cols));
    for (int ii = 0; ii < N; ++ii) {
      for (int i = 0; i < M2; ++i) {
        for (int j = 0; j < experiment->inputDim; ++j) {
          double entry = X_ML[ii](i, j);
          file.write(reinterpret_cast<char *>(&entry), sizeof(double));
        }
      }
    }
    file.close();
  }
}

void ExpectedInformationEstimator::Allocate()
{
  cout << "ExpectedInformationEstimator::Allocate" << endl;

  X_prior = MatrixXd::Zero(N, experiment->inputDim);
  G_prior = MatrixXd::Zero(N, experiment->outputDim);
  Y_prior = MatrixXd::Zero(N, experiment->outputDim);

  priorLogDensities   = VectorXd::Zero(N);
  priorLogLikelihoods = MatrixXd::Zero(N, N);

  q_ML.resize(N);
  q_CL.resize(N);

  postMean.resize(N);
  postCov.resize(N);
  condMean.resize(N);
  condCov.resize(N);

  X_ML.resize(N, MatrixXd::Zero(M1, experiment->inputDim));
  G_ML.resize(N, MatrixXd::Zero(M1, experiment->outputDim));
  X_CL.resize(N, MatrixXd::Zero(M1, experiment->inputDim));
  G_CL.resize(N, MatrixXd::Zero(M1, experiment->outputDim));

  logf_ML    = MatrixXd::Zero(N, M1);
  logw_ML    = MatrixXd::Zero(N, M1);
  logML      = VectorXd::Zero(N);
  logMLexact = VectorXd::Zero(N);

  logf_CL    = MatrixXd::Zero(N, M2);
  logw_CL    = MatrixXd::Zero(N, M2);
  logCL      = VectorXd::Zero(N);
  logCLexact = VectorXd::Zero(N);

  ESS = VectorXd::Zero(N);

  if (useMIS) {
    AllocateMIS();
    mixtureIndices.resize(N);
  }

  allocated = true;
}

void ExpectedInformationEstimator::AllocateMIS()
{
  // int Lmax = N + (N - 1) * M1;
  ML_priorLogDensities = MatrixXd::Zero(N, M1);
  ML_biasLogDensities  = MatrixXd::Zero(N, M1);
  allocatedMIS         = true;
}

// Sort the samples of X_prior before evaluating in descending order of prior density
void ExpectedInformationEstimator::SortPriorSamples()
{
  cout << "ExpectedInformationEstimator::SortPriorSamples" << endl;
  vector<tuple<double, int>> pairs;
  for (int i = 0; i < N; ++i) {
    pairs.emplace_back(priorLogDensities(i), i);
  }
  sort(pairs.rbegin(), pairs.rend());
  MatrixXd original          = X_prior;
  VectorXd originalDensities = priorLogDensities;
  for (int i = 0; i < N; ++i) {
    X_prior.row(i)       = original.row(get<1>(pairs[i]));
    priorLogDensities(i) = originalDensities(get<1>(pairs[i]));
  }
}

void ExpectedInformationEstimator::EvaluatePriorSamples(const VectorXd& design)
{
  // Evaluate the prior outputs and observations
  for (int i = 0; i < N; ++i) {
    G_prior.row(i) = experiment->Evaluate(X_prior.row(i), design);
    Y_prior.row(i) = G_prior.row(i) + experiment->GetNoiseSample().transpose();
  }
}

double ExpectedInformationEstimator::Evaluate(const VectorXd& design)
{
  cout << "ExpectedInformationEstimator::Evaluate" << endl;
  Timer timer;
  timer.Start();

  if (allocated == false) {
    Allocate();
    if (useMIS) {
      AllocateMIS();
    }
  }

  if (reusePriorSamples == false) {
    for (int i = 0; i < N; ++i) {
      X_prior.row(i)       = experiment->GetPriorSample();
      priorLogDensities(i) = experiment->PriorLogDensity(X_prior.row(i));
    }
    if (useSortHeuristic) {
      SortPriorSamples();
    }
    EvaluatePriorSamples(design);
  }

  double EIG = 0;

  if (useMIS) {
    EIG = EvaluateMIS(design);
  } else if (useIS)   {
    EIG = EvaluateIS(design);
  } else {
    EIG = EvaluateNoIS(design);
  }

  executionTime = timer.Elapsed();

  return EIG;
}

double ExpectedInformationEstimator::EvaluateBiasCorrectedEIG()
{
  double EIG = 0;
  for (int i = 0; i < N; ++i) {
    VectorXd CL         = exp((logf_CL.row(i) + logw_CL.row(i)).array());
    VectorXd ML         = exp((logf_ML.row(i) + logw_ML.row(i)).array());
    double   mu_CL      = exp(logCL(i));
    double   mu_ML      = exp(logML(i));
    double   var_CL     = (CL - VectorXd::Constant(M2, mu_CL)).squaredNorm() / (M2 - 1);
    double   var_ML     = (ML - VectorXd::Constant(M1, mu_ML)).squaredNorm() / (M1 - 1);
    double   correction = 0.5 * (var_CL / (mu_CL * mu_CL * M2) - var_ML / (mu_ML * mu_ML * M1));
    EIG += logCL(i) - logML(i) + correction;
  }
  return EIG / N;
}

double ExpectedInformationEstimator::EvaluateNoIS(const VectorXd& design)
{
  cout << "ExpectedInformationEstimator::EvaluateNoIS" << endl;
  // Outer Monte Carlo
  for (int i = 0; i < N; ++i) {
    // Evaluate log marginal likelihood
    for (int j = 0; j < M1; ++j) {
      X_ML[i].row(j) = experiment->GetPriorSample();
      G_ML[i].row(j) = experiment->Evaluate(X_ML[i].row(j), design);
      logf_ML(i, j)  = experiment->LogLikelihoodCached(Y_prior.row(i), X_ML[i].row(j), G_ML[i].row(j), design);
    }
    logML(i) = LogSumExp(logf_ML.row(i)) - log(M1);

    // Evaluate log conditional likelihood
    if (experiment->nuisanceParameters == 0) {
      logCL(i) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(i), G_prior.row(i), design);
    } else {
      for (int k = 0; k < M2; ++k) {
        X_CL[i].row(k)                                      = X_prior.row(i);
        X_CL[i].row(k).tail(experiment->nuisanceParameters) = experiment->GetCondPriorSample(X_prior.row(i));
        G_CL[i].row(k)                                      = experiment->Evaluate(X_CL[i].row(k), design);
        logf_CL(i, k)                                       = experiment->LogLikelihoodCached(Y_prior.row(i), X_CL[i].row(k), G_CL[i].row(k), design);
      }
      logCL(i) = LogSumExp(logf_CL.row(i)) - log(M2);
    }
  }

  if (useBiasCorrection) {
    return EvaluateBiasCorrectedEIG();
  } else {
    return (logCL - logML).array().sum() / N;
  }
}

double ExpectedInformationEstimator::EvaluateIS(const VectorXd& design)
{
  cout << "ExpectedInformationEstimator::EvaluateIS" << endl;
  // evaluate log likelihood matrix with (i,j) entries log p(y^(i)|X_prior(j),d)
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      priorLogLikelihoods(i, j) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(j), G_prior.row(j), design);
    }
    priorLogLikelihoods.row(i).array() -= LogSumExp(priorLogLikelihoods.row(i)); // normalize
  }

  // outer Monte Carlo index
  for (int i = 0; i < N; ++i) {
    if (useExactPosterior) {
      postMean[i] = experiment->PosteriorMean(Y_prior.row(i), design);
      postCov[i]  = experiment->PosteriorCovariance(Y_prior.row(i), design);
    } else {
      VectorXd importanceWeights = priorLogLikelihoods.row(i).array().exp();
      ESS(i) = EffectiveSampleSize(importanceWeights);

      // compute posterior mean and covariance as weighted sample mean and covariance
      postMean[i] = WeightedSampleMean(X_prior, importanceWeights);
      postCov[i]  = WeightedSampleCovariance(X_prior, importanceWeights);
      postCov[i] += nugget * MatrixXd::Identity(experiment->inputDim, experiment->inputDim);
    }

    // compute conditional mean and covariance assuming MVN
    if (experiment->nuisanceParameters > 0) {
      tie(condMean[i], condCov[i]) = MultivariateNormal::GetConditional(postMean[i], postCov[i], X_prior.row(i).head(experiment->parametersOfInterest));
      // condCov[i]               += nugget * MatrixXd::Identity(experiment->nuisanceParameters, experiment->nuisanceParameters);
    }

    // construct biasing distributions q_ML and q_CL as MVN or MVT
    if (biasingDistributionType == "MVN") {
      q_ML[i] = make_shared<MultivariateNormal>(postMean[i], postCov[i]);
      if (experiment->nuisanceParameters > 0) {
        q_CL[i] = make_shared<MultivariateNormal>(condMean[i], condCov[i]);
      }
    } else if (biasingDistributionType == "MVT")   {
      q_ML[i] = make_shared<MultivariateT>(postMean[i], postCov[i], dof);
      if (experiment->nuisanceParameters > 0) {
        q_CL[i] = make_shared<MultivariateT>(condMean[i], condCov[i], dof);
      }
    }

    // evaluate log marginal likelihood estimator
    for (int j = 0; j < M1; ++j) {
      X_ML[i].row(j) = q_ML[i]->GetSample();
      G_ML[i].row(j) = experiment->Evaluate(X_ML[i].row(j), design);
      logf_ML(i, j)  = experiment->LogLikelihoodCached(Y_prior.row(i), X_ML[i].row(j), G_ML[i].row(j), design);
      logw_ML(i, j)  = experiment->PriorLogDensity(X_ML[i].row(j)) - q_ML[i]->LogDensity(X_ML[i].row(j));
    }
    logML(i) = LogSumExp(logf_ML.row(i) + logw_ML.row(i)) - log(M1);

    // evaluate log conditional likelihood estimator
    if (experiment->nuisanceParameters == 0) {
      logCL(i) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(i), G_prior.row(i), design);
    } else {
      for (int k = 0; k < M2; ++k) {
        X_CL[i].row(k)                                      = X_prior.row(i);
        X_CL[i].row(k).tail(experiment->nuisanceParameters) = q_CL[i]->GetSample();
        G_CL[i].row(k)                                      = experiment->Evaluate(X_CL[i].row(k), design);
        logf_CL(i, k)                                       = experiment->LogLikelihoodCached(Y_prior.row(i), X_CL[i].row(k), G_CL[i].row(k), design);
        logw_CL(i, k)                                       = experiment->CondPriorLogDensity(X_CL[i].row(k)) - q_CL[i]->LogDensity(X_CL[i].row(k).tail(experiment->nuisanceParameters));
      }
      logCL(i) = LogSumExp(logf_CL.row(i) + logw_CL.row(i)) - log(M2);
    }
  }
  if (useBiasCorrection) {
    return EvaluateBiasCorrectedEIG();
  } else {
    return (logCL - logML).array().sum() / N;
  }
}

double ExpectedInformationEstimator::EvaluateMIS(const VectorXd& design)
{
  cout << "ExpectedInformationEstimator::EvaluateMIS" << endl;
  // evaluate log likelihood matrix with (i,j) entries log p(y^(i)|X_prior(j),d)
  // all of these will be used when computing the mixtureLogLikelihoods
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      priorLogLikelihoods(i, j) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(j), G_prior.row(j), design);
    }
    priorLogLikelihoods.row(i).array() -= LogSumExp(priorLogLikelihoods.row(i)); // normalize
  }

  int rescaled = 0;

  // outer Monte Carlo index
  for (int i = 0; i < N; ++i) {
    // initialize mixture with prior (index -1, just in case we try to access it)
    mixtureIndices[i].resize(1);
    mixtureIndices[i][0] = -1;

    // vector<int> mixtureIndices(1, -1);
    int mixtureSize = 1;

    // identify most useful q_ML from previous iterations to include in the mixture
    if (i > 0) {
      if (useReverseLikelihood) {
        // sort log reverse likelihoods in descending order
        vector<tuple<double, int>> pairs;
        for (int m = 0; m < i; ++m) {
          pairs.emplace_back(q_ML[m]->LogDensity(X_prior.row(i)), m);
        }
        sort(pairs.rbegin(), pairs.rend());

        // only include q_ML where reverse likelihood > prior density
        for (int m = 0; m < min(i, maxComponents - 1); ++m) {
          int index = get<1>(pairs[m]);
          if (get<0>(pairs[m]) > priorLogDensities(index)) {
            mixtureIndices[i].push_back(index);
          }
        }
      }
      if (useMinSampleDistance) {
        // sort min distance in ascending order
        vector<tuple<double, int>> pairs;
        for (int m = 0; m < i; ++m) {
          double dist, minDist = (X_ML[m].row(0) - X_prior.row(i)).squaredNorm();
          for (int j = 1; j < M1; j++) {
            dist = (X_ML[m].row(j) - X_prior.row(i)).squaredNorm();
            if (dist < minDist) {
              minDist = dist;
            }
          }
          pairs.emplace_back(minDist, m);
        }
        sort(pairs.begin(), pairs.end());
        for (int m = 0; m < min(i, maxComponents - 1); ++m) {
          int index = get<1>(pairs[m]);
          mixtureIndices[i].push_back(index);
        }
      }

      mixtureSize = mixtureIndices[i].size();
    }

    // number of mixture samples in this iteration
    int L = N + (mixtureSize - 1) * M1;

    // populate rows of mixtureSamples with L mixture samples
    // where the first N rows are from the prior
    // and every subsequent block of M1 rows are from q_ML[mixtureIndices[i][m]]
    MatrixXd mixtureSamples(L, experiment->inputDim);
    mixtureSamples.topRows(N) = X_prior;
    for (int m = 1; m < mixtureSize; ++m) {
      mixtureSamples.block(N + (m - 1) * M1, 0, M1, experiment->inputDim) = X_ML[mixtureIndices[i][m]];
    }

    // populate columns of mixtureLogDensities with mixture components evaluated at the L mixture samples
    MatrixXd mixtureLogDensities(L, mixtureSize);

    // first segment of N elements of first column were previously calculated in Evaluate()
    mixtureLogDensities.block(0, 0, N, 1) = priorLogDensities;

    // remaining mixtureSize-1 segments of M1 elements in the first column were
    // previously calculated as part of the importance weights when estimating the marginal likelihood
    for (int m = 1; m < mixtureSize; ++m) {
      mixtureLogDensities.block(N + (m - 1) * M1, 0, M1, 1) = ML_priorLogDensities.row(mixtureIndices[i][m]).transpose();
    }

    // populate the remaining columns
    for (int m = 1; m < mixtureSize; ++m) {
      // the first N elements of each column are q_ML[mixtureIndicies[m]] evaluated at samples from X_prior
      // check if q_ML[mixtureIndices[i][m]] has been evaluated at the prior samples
      if (ML_priorCache.count(mixtureIndices[i][m]) == 0) {
        // evaluate q_ML[mixtureIndices[m]] at each of the N prior samples in X_prior
        ML_priorCache[mixtureIndices[i][m]] = VectorXd::Zero(N);
        for (int ii = 0; ii < N; ++ii) {
          ML_priorCache[mixtureIndices[i][m]](ii) = q_ML[mixtureIndices[i][m]]->LogDensity(X_prior.row(ii));
        }
      }
      mixtureLogDensities.block(0, m, N, 1) = ML_priorCache[mixtureIndices[i][m]];

      // the remaining mixtureSize-1 segments of length M1 in column m correspond to
      // q_ML[mixtureIndices[i][m]] evaluated at the M1 samples in X_ML[mixtureIndices[i][n]]
      for (int n = 1; n < mixtureSize; ++n) {
        // tuple(a,b) -> q_ML[a] evaluated at X_ML[b]
        tuple<int, int> t(mixtureIndices[i][m], mixtureIndices[i][n]);

        // q_ML[mixtureIndices[i][m]] has already been evaluated at X_ML[[mixtureIndices[m]]
        // when computing the marginal likelihood in previous iterations
        if (n == m) {
          ML_cache[t] = ML_biasLogDensities.row(mixtureIndices[i][m]).transpose();
        }

        // check to see if the tuple corresponds to an existing entry in the cache
        if (ML_cache.count(t) == 0) {
          // evaluate q_ML[mixtureIndices[i][m]] at the M1 samples in X_ML[mixtureIndices[i][n]]
          ML_cache[t] = VectorXd::Zero(M1);
          for (int j = 0; j < M1; ++j) {
            ML_cache[t](j) = q_ML[mixtureIndices[i][m]]->LogDensity(X_ML[mixtureIndices[i][n]].row(j));
          }
        }
        mixtureLogDensities.block(N + (n - 1) * M1, m, M1, 1) = ML_cache[t];
      }
    }

    // likelihood of y(i) given the L samples in the mixture
    VectorXd mixtureLogLikelihoods(L);

    // the first N elements have been previously computed in the NxN priorLogLikelihoods matrix
    mixtureLogLikelihoods.head(N) = priorLogLikelihoods.row(i).transpose();

    // the remaning mixtureSize-1 segments of length M1 have not been evaluated - no cache
    for (int m = 1; m < mixtureSize; ++m) {
      // evaluate the log likelihood of y(i) at the M1 samples of X_ML[mixtureIndices[i][m]] using the
      // cached model evaluations G_ML[mixtureIndices[i][m]]
      for (int j = 0; j < M1; ++j) {
        mixtureLogLikelihoods(N + (m - 1) * M1 + j) = experiment->LogLikelihoodCached(Y_prior.row(i), X_ML[mixtureIndices[i][m]].row(j), G_ML[mixtureIndices[i][m]].row(j), design);
      }
    }

    // mixture weights, N/L for the first component (prior), and M1/L for remaining components
    VectorXd logMixtureWeights(mixtureSize);
    logMixtureWeights    = VectorXd::Ones(mixtureSize) * (log(M1) - log(L));
    logMixtureWeights[0] = log(N) - log(L);

    // MIS importance weights computed as likelihood * prior / q_MIS
    // note: column 0 of mixtureLogDensities corresponds to first component (prior)
    VectorXd logImportanceWeights(L);
    for (int ell = 0; ell < L; ++ell) {
      logImportanceWeights(ell) = mixtureLogLikelihoods(ell) + mixtureLogDensities(ell, 0) - LogSumExp(mixtureLogDensities.row(ell).transpose() + logMixtureWeights);
    }

    // normalize the log weights by subtracting log normalizing constant
    VectorXd importanceWeights = (logImportanceWeights.array() - LogSumExp(logImportanceWeights)).exp();

    // compute effective sample size
    ESS(i) = EffectiveSampleSize(importanceWeights);

    // compute posterior mean and covariance as weighted sample mean and covariance
    postMean[i] = WeightedSampleMean(mixtureSamples, importanceWeights);
    postCov[i]  = WeightedSampleCovariance(mixtureSamples, importanceWeights);
    postCov[i] += nugget * MatrixXd::Identity(experiment->inputDim, experiment->inputDim);

    // compute conditional mean and covariance assuming MVN
    if (experiment->nuisanceParameters > 0) {
      if (useMarginal) {
        condMean[i] = postMean[i].tail(experiment->nuisanceParameters);
        condCov[i]  = postCov[i].bottomRightCorner(experiment->nuisanceParameters, experiment->nuisanceParameters);
      } else {
        tie(condMean[i], condCov[i]) = MultivariateNormal::GetConditional(postMean[i], postCov[i], X_prior.row(i).head(experiment->parametersOfInterest));
      }
    }

    // construct biasing distributions q_ML^(i) and q_CL^(i)
    if (biasingDistributionType == "MVN") {
      q_ML[i] = make_shared<MultivariateNormal>(postMean[i], postCov[i]);
      if (experiment->nuisanceParameters > 0) {
        q_CL[i] = make_shared<MultivariateNormal>(condMean[i], condCov[i]);
      }
    } else if (biasingDistributionType == "MVT")   {
      q_ML[i] = make_shared<MultivariateT>(postMean[i], postCov[i], dof);
      if (experiment->nuisanceParameters > 0) {
        q_CL[i] = make_shared<MultivariateT>(condMean[i], condCov[i], dof);
      }
    }

    // evaluate log marginal likelihood estimator
    for (int j = 0; j < M1; ++j) {
      X_ML[i].row(j)             = q_ML[i]->GetSample();
      G_ML[i].row(j)             = experiment->Evaluate(X_ML[i].row(j), design);
      ML_priorLogDensities(i, j) = experiment->PriorLogDensity(X_ML[i].row(j));
      ML_biasLogDensities(i, j)  = q_ML[i]->LogDensity(X_ML[i].row(j));
      logf_ML(i, j)              = experiment->LogLikelihoodCached(Y_prior.row(i), X_ML[i].row(j), G_ML[i].row(j), design);
      logw_ML(i, j)              = ML_priorLogDensities(i, j) - ML_biasLogDensities(i, j);
    }
    logML(i) = LogSumExp(logf_ML.row(i) + logw_ML.row(i)) - log(M1);

    if (debugCondLikelihood) {
      // Evaluate log conditional likelihood without using IS
      if (experiment->nuisanceParameters == 0) {
        logCL(i) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(i), G_prior.row(i), design);
      } else {
        for (int k = 0; k < M2; ++k) {
          X_CL[i].row(k)                                      = X_prior.row(i);
          X_CL[i].row(k).tail(experiment->nuisanceParameters) = experiment->GetCondPriorSample(X_prior.row(i));
          G_CL[i].row(k)                                      = experiment->Evaluate(X_CL[i].row(k), design);
          logf_CL(i, k)                                       = experiment->LogLikelihoodCached(Y_prior.row(i), X_CL[i].row(k), G_CL[i].row(k), design);
        }
        logCL(i) = LogSumExp(logf_CL.row(i)) - log(M2);
      }
    } else {
      // evaluate log conditional likelihood estimator
      if (experiment->nuisanceParameters == 0) {
        logCL(i) = experiment->LogLikelihoodCached(Y_prior.row(i), X_prior.row(i), G_prior.row(i), design);
      } else {
        for (int k = 0; k < M2; ++k) {
          X_CL[i].row(k)                                      = X_prior.row(i);
          X_CL[i].row(k).tail(experiment->nuisanceParameters) = q_CL[i]->GetSample();
          G_CL[i].row(k)                                      = experiment->Evaluate(X_CL[i].row(k), design);
          logf_CL(i, k)                                       = experiment->LogLikelihoodCached(Y_prior.row(i), X_CL[i].row(k), G_CL[i].row(k), design);
          logw_CL(i, k)                                       = experiment->CondPriorLogDensity(X_CL[i].row(k)) - q_CL[i]->LogDensity(X_CL[i].row(k).tail(experiment->nuisanceParameters));
        }
        logCL(i) = LogSumExp(logf_CL.row(i) + logw_CL.row(i)) - log(M2);
      }
    }
  }

  double EIG = useBiasCorrection ? EvaluateBiasCorrectedEIG() : (logCL - logML).array().sum() / N;

  if (std::isnan(EIG)) {
    cout << "encountered NaN" << endl;
    WriteToFile("failure", design);
  }

  return EIG;
}
