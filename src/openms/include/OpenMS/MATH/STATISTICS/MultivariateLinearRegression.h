// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/RegressionUtils.h>

#include <Eigen/Dense>

using namespace Eigen;

namespace OpenMS
{
  namespace Math
  {
    /*
      @brief Estimates the parameters of a linear equation mapping multiple variables to a single value.

      The equation will be of the form:
      y = a + b*x + c*t + d*s ...
    */
    class OPENMS_DLLAPI MultivariateLinearRegression
    {
public:
      MultivariateLinearRegression();

      /**
          @brief Compute a linear regression using multiple variables
      */
      void computeRegression(std::vector<double> target, std::list<std::vector<double>> variables);

      /// Access to the parameters computed in the regression
      std::vector<double> getRegressionParams() const;

      /// Relative error of the fit
      double getRelativeError() const;

      std::vector<double> getPredictedValues();

protected:
        std::vector<double> regression_params_;
        double relative_error_;
        std::vector<double> predicted_values_;
    }; //class

    void MultivariateLinearRegression::computeRegression(std::vector<double> target, std::list<std::vector<double>> variables)
    {

      int num_rows = target.size();
      // create a vector of ones to put the matrix to calculate the intercept
      VectorXd ones = VectorXd::Ones(num_rows);

      MatrixXd A(num_rows, variables.size());
      A << ones;
      for (auto var : variables)
      {
          A << Eigen::Map<VectorXd>(var.data(), num_rows);
      }
      /**Eigen::Map<VectorXd>(theo_im.data(), num_rows),
      Eigen::Map<VectorXd>(mz.data(), num_rows),
      Eigen::Map<VectorXd>(charges.data(), num_rows),
      Eigen::Map<VectorXd>(rt.data(), num_rows);
    */
      VectorXd b = Eigen::Map<VectorXd>(target.data(), num_rows);
      // this should return a vector of parameters by computing the singular value decomposition
      Eigen::VectorXd regression_params = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
      Eigen::VectorXd predicted_values = A*regression_params;

      regression_params_ = std::vector<double>(regression_params.data(), regression_params.data()+regression_params.size());
      predicted_values_ = std::vector<double>(predicted_values.data(), predicted_values.data()+predicted_values.size());
      relative_error_ = (predicted_values - b).norm() / b.norm();

    }

  } // namespace Math
} // namespace OpenMS


