#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <Eigen/Core>
#include "Utilities.h"
#include "MultivariateNormal.h"
#include "RandomGenerator.h"
#include "merr/PolynomialModel.h"
#include "pce/PCExpansion.h"
#include "MCMC.h"

using namespace std;
using namespace Eigen;
using namespace Utilities;
using namespace merr;

int  RunTests(int argc, char **argv);
void RunPosteriorGrid(shared_ptr<Model> model, const string& xgridfile, const string& pgridfile);
void RunMCMC(shared_ptr<Model> model, const int nmcmc, const string& mcmcfile);

int  main(int argc, char **argv)
{
  RunTests(argc, argv);

  cout << "\nModelInference\n" << endl;

  string xdatafile      = GetOption<string>(argc, argv, "-xdatafile", "none");
  string ydatafile      = GetOption<string>(argc, argv, "-ydatafile", "none");
  string modeltype      = GetOption<string>(argc, argv, "-modeltype", "constant");
  string likelihoodtype = GetOption<string>(argc, argv, "-likelihoodtype", "full");
  int    nsamps         = GetOption<int>(argc, argv, "-nsamps", 1000);
  bool   usenisp        = GetOption<bool>(argc, argv, "-usenisp", false);
  int    pcdim          = GetOption<int>(argc, argv, "-pcdim", 5);
  double abceps         = GetOption<double>(argc, argv, "-abceps", 0.1);
  double abcgamma       = GetOption<double>(argc, argv, "-abceps", 1.0);
  string priortype      = GetOption<string>(argc, argv, "-priortype", "uniform");
  double datanoise      = GetOption<double>(argc, argv, "-datanoise", 0.1);
  bool   infernoise     = GetOption<bool>(argc, argv, "-infernoise", false);
  string xgridfile      = GetOption<string>(argc, argv, "-xgridfile", "none");
  string paramsfile     = GetOption<string>(argc, argv, "-paramsfile", "none");
  string pgridfile      = GetOption<string>(argc, argv, "-pgridfile", "none");
  string mcmcfile       = GetOption<string>(argc, argv, "-mcmcfile", "none");
  int    nmcmc          = GetOption<int>(argc, argv, "-nmcmc", 0);

  if ((xdatafile == "none") || (ydatafile == "none")) {
    cout << "No data files, exiting." << endl;
    return 0;
  }

  shared_ptr<Model> model;

  if (modeltype == "constant") model = make_shared<PolynomialModel<0>>();
  if (modeltype == "linear") model = make_shared<PolynomialModel<1>>();
  if (modeltype == "quadratic") model = make_shared<PolynomialModel<2>>();
  if (modeltype == "cubic") model = make_shared<PolynomialModel<3>>();

  if (likelihoodtype == "abc") model->likelihood = Model::LikelihoodType::ABC;
  if (likelihoodtype == "gauss") model->likelihood = Model::LikelihoodType::Gauss;
  if (likelihoodtype == "marginal") model->likelihood = Model::LikelihoodType::Marginal;
  if (likelihoodtype == "gaussmarginal") model->likelihood = Model::LikelihoodType::GaussMarginal;
  if (likelihoodtype == "full") model->likelihood = Model::LikelihoodType::Full;

  model->name       = modeltype;
  model->inferNoise = infernoise;
  model->nsamps     = nsamps;
  model->useNISP    = usenisp;
  model->pcdim      = pcdim;
  model->abceps     = abceps;
  model->abcgamma   = abcgamma;
  model->datanoise  = datanoise;
  model->SetData(ReadEigenAsciiFile(xdatafile), ReadEigenAsciiFile(ydatafile));

  cout << model << endl;

  if ((xgridfile != "none") && (pgridfile != "none")) {
    RunPosteriorGrid(model, xgridfile, pgridfile);
  }

  if ((mcmcfile != "none") && (nmcmc > 0)) {
    RunMCMC(model, nmcmc, mcmcfile);
  }

  return 0;
}

void RunMCMC(shared_ptr<Model> model, const int nmcmc, const string& mcmcfile)
{
  cout << "Evaluating MCMC" << endl;
  auto f = [&](const VectorXd& alpha) {
             return model->LogPosterior(alpha);
           };
  auto mcmc = make_shared<MCMC>(f, model->GetMCMCDim());
  mcmc->SetAdaptBeta(0.1);
  mcmc->Run(nmcmc);
  auto chain = mcmc->GetChain();
  cout << "post mean" << chain.colwise().mean() << endl;
  cout << "Writing MCMC results to " << mcmcfile << ".*" << endl;
  WriteEigenAsciiFile(mcmcfile + ".chain", mcmc->GetChain());
  WriteEigenAsciiFile(mcmcfile + ".logdens", mcmc->GetLogDensities());
  WriteEigenAsciiFile(mcmcfile + ".acf", mcmc->GetAutocorrelation(100));
}

void RunPosteriorGrid(shared_ptr<Model> model, const string& xgridfile, const string& pgridfile)
{
  cout << "Evaluating posterior on grid" << endl;
  MatrixXd inputs        = ReadEigenAsciiFile(xgridfile);
  VectorXd logPosteriors = VectorXd::Zero(inputs.rows());
  // #pragma omp parallel for schedule (dynamic)
  for (int i = 0; i < inputs.rows(); ++i) {
    if ((i + 1) % (inputs.rows() / 20) == 0) {
      cout << "."; cout.flush();
    }
    logPosteriors(i) = model->LogPosterior(inputs.row(i));
  }
  cout << "\nWriting posterior values to " << pgridfile << endl;
  ofstream out(pgridfile);
  for (int i = 0; i < inputs.rows(); ++i) {
    out << inputs.row(i) << " " << logPosteriors(i) << endl;
  }
  out.close();
}

int RunTests(int argc, char **argv)
{
  auto f = [](const VectorXd& x) {
             return x(0) * x(1) + pow(x(1), 2) + pow(x(0), 3);
           };
  auto pce = make_shared<PCExpansion>();
  pce->SetFunction(f);
  pce->SetVariables({ PCExpansion::GermType::Normal, PCExpansion::GermType::Normal });
  pce->SetMaxOrder(4);
  auto moments = pce->GetMoments();

  cout << "mean " << moments.first << endl;
  cout << "var  " << moments.second << endl;

  double mse = 0;
  for (int i = 0; i < 1000; ++i) {
    VectorXd xi = RandomGenerator::GetNormalRandomVector(2);
    mse += pow(f(xi) - pce->Evaluate(xi), 2);
  }
  cout << "mse " << (mse / 1000) << endl;

  int  dim            = 2;
  MatrixXd covariance = MatrixXd::Identity(dim, dim) * 0.2;
  auto mvn1           = make_shared<MultivariateNormal>(VectorXd::Ones(dim), covariance);
  auto mvn2           = make_shared<MultivariateNormal>(VectorXd::Zero(dim), covariance);
  auto logdens        = [&](const VectorXd& sample) {
                          return log(exp(mvn1->LogDensity(sample)) + exp(mvn2->LogDensity(sample)));
                        };
  auto mcmc = make_shared<MCMC>(logdens, dim);
  mcmc->Run(10000);
  WriteEigenAsciiFile("test.chain", mcmc->GetChain());
  WriteEigenAsciiFile("test.logdens", mcmc->GetLogDensities());
  WriteEigenAsciiFile("test.acf", mcmc->GetAutocorrelation(100));

  return 0;
}

