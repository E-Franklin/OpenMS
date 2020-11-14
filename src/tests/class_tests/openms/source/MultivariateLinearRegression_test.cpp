// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/MultivariateLinearRegression.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace Math;
using namespace Eigen;

START_TEST(MultivariateLinearRegression, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MultivariateLinearRegression* ptr = nullptr;
MultivariateLinearRegression* null_ptr = nullptr;
START_SECTION(MultivariateLinearRegression())
{
	ptr = new MultivariateLinearRegression();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MultivariateLinearRegression())
{
	delete ptr;
}
END_SECTION

// Create a test data input for the equation Ax = b
// The columns for A are passed as vectors, x is the expected params and b is the target
// Test with a single variable
vector<double> x1(10);
vector<double> x2 {1,2,7,4,5,6,7,8,9,10};
vector<double> target(10);
vector<double> target1(10);
vector<double> target2(10);
vector<double> multivar_target(10);
vector<double> multivar_target1(10);
vector<double> expected_params({4, 5.5});
vector<double> expected_params1({0, 5.5});
vector<double> expected_params2({4.81818, 5.41818});
vector<double> multivar_params({4, 5.5, 7});
vector<double> multivar_params1({1.725, 5.325, 7.375});
vector<double> expected_predicted(10);
vector<double> multivar_expected_predicted(10);

for (int i=0; i < 10; ++i)
{
  x1[i] = i;

  target[i] = expected_params[0] + expected_params[1]*x1[i];
  target1[i] = expected_params1[1]*x1[i]; // no intercept
  target2[i] = expected_params[0] + expected_params[1]*x1[i];

  expected_predicted[i] = expected_params2[1]*i + expected_params2[0];

  multivar_target[i] = 4 + 5.5*x1[i] + 7*x2[i];
  multivar_target1[i] = 4 + 5.5*x1[i] + 7*x2[i];

  multivar_expected_predicted[i] = multivar_params1[0] + multivar_params1[1] * x1[i] + multivar_params1[2] * x2[i];
}

target2[3] = 25;
multivar_target1[3] = 40;

MultivariateLinearRegression m_reg;

START_SECTION((void computeRegression(std::vector<double> target, std::vector<std::vector<double>> variables)))
{
  // test with a single variable, with and without intercept, with no error
  m_reg.computeRegression(target, {x1});
  vector<double> calculated_params = m_reg.getRegressionParams();
  vector<double> predicted = m_reg.getPredictedValues();
  TEST_REAL_SIMILAR(calculated_params[0], expected_params[0])
  TEST_REAL_SIMILAR(calculated_params[1], expected_params[1])
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.0)
  for (int i=0; i < 10; ++i)
  {
    TEST_REAL_SIMILAR(predicted[i], target[i])
  }

  m_reg.computeRegression(target1, {x1});
  calculated_params = m_reg.getRegressionParams();
  predicted = m_reg.getPredictedValues();
  TEST_REAL_SIMILAR(calculated_params[0], expected_params1[0])
  TEST_REAL_SIMILAR(calculated_params[1], expected_params1[1])
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.0)
  for (int i=0; i < 10; ++i)
  {
    TEST_REAL_SIMILAR(predicted[i], target1[i])
  }

  // test with a single variable, with intercept, and with error
  m_reg.computeRegression(target2, {x1});
  calculated_params = m_reg.getRegressionParams();
  predicted = m_reg.getPredictedValues();
  TEST_REAL_SIMILAR(calculated_params[0], expected_params2[0])
  TEST_REAL_SIMILAR(calculated_params[1], expected_params2[1])
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.040144)
  for (int i=0; i < 10; ++i)
  {
    TEST_REAL_SIMILAR(predicted[i], expected_predicted[i])
  }

  // test with three variables, with intercept, with and without error
  m_reg.computeRegression(multivar_target, {x1, x2});
  calculated_params = m_reg.getRegressionParams();
  predicted = m_reg.getPredictedValues();
  TEST_REAL_SIMILAR(calculated_params[0], multivar_params[0])
  TEST_REAL_SIMILAR(calculated_params[1], multivar_params[1])
  TEST_REAL_SIMILAR(calculated_params[2], multivar_params[2])
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.0)
  for (int i=0; i < 10; ++i)
  {
    TEST_REAL_SIMILAR(predicted[i], multivar_target[i])
  }

  m_reg.computeRegression(multivar_target1, {x1, x2});
  calculated_params = m_reg.getRegressionParams();
  predicted = m_reg.getPredictedValues();
  TEST_REAL_SIMILAR(calculated_params[0], multivar_params1[0])
  TEST_REAL_SIMILAR(calculated_params[1], multivar_params1[1])
  TEST_REAL_SIMILAR(calculated_params[2], multivar_params1[2])
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.031894)
  for (int i=0; i < 10; ++i)
  {
    TEST_REAL_SIMILAR(predicted[i], multivar_expected_predicted[i])
  }
}
END_SECTION

START_SECTION((double getRegressionParams() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getRelativeError() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getPredictedValues() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



