// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{
  namespace MapConversion
  {
    void convert(UInt64 const input_map_index,
                 MSExperiment<> & input_map,
                 ConsensusMap & output_map,
                 Size n)
    {
      output_map.clear(true);

      // see @todo above
      output_map.setUniqueId();

      input_map.updateRanges(1);
      if (n > input_map.getSize())
      {
        n = input_map.getSize();
      }
      output_map.reserve(n);
      std::vector<Peak2D> tmp;
      tmp.reserve(input_map.getSize());

      // TODO Avoid tripling the memory consumption by this call
      input_map.get2DData(tmp);

      std::partial_sort(tmp.begin(),
                        tmp.begin() + n,
                        tmp.end(),
                        reverseComparator(Peak2D::IntensityLess()));

      for (Size element_index = 0; element_index < n; ++element_index)
      {
        output_map.push_back(ConsensusFeature(input_map_index,
                                              tmp[element_index],
                                              element_index));
      }

      output_map.getFileDescriptions()[input_map_index].size = n;
      output_map.updateRanges();
    }
  } // namespace MapConversion
} // namespace OpenMS

