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
vector<double> var(10);
vector<double> target(10);
vector<double> target1(10);
vector<double> expected_params({4, 5.5});
vector<double> expected_params1({0, 5.5});

for (int i=0; i < 10; ++i)
{
  var[i] = i;
  target[i] = 5.5*i + 4;
  target1[i] = 5.5*i; // no intercept
}

MultivariateLinearRegression m_reg;

START_SECTION((void computeRegression(std::vector<double> target, std::vector<std::vector<double>> variables)))
{
  // test with a single variable, with and without intercept, with no error
  m_reg.computeRegression(target, {var});
  TEST_EQUAL(m_reg.getRegressionParams()==expected_params, true)
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.0)
  TEST_EQUAL(m_reg.getPredictedValues()==target, true)

  m_reg.computeRegression(target1, {var});
  TEST_EQUAL(m_reg.getRegressionParams()==expected_params1, true)
  TEST_REAL_SIMILAR(m_reg.getRelativeError(), 0.0)
  TEST_EQUAL(m_reg.getPredictedValues()==target1, true)
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



