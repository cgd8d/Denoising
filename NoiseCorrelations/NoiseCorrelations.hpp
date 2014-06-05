#ifndef NoiseCorrelations_hpp
#define NoiseCorrelations_hpp

#include "Utilities/IndexHandler.hpp"
#include "EXOUtilities/SystemOfUnits.hh"
#include <vector>

/*
Class for holding noise correlation information at a particular frequency.
We store this information in memory as a single unpacked matrix;
the matrix is guaranteed to be symmetric.

We store real and imaginary correlations always,
even though the imaginary-imaginary and imaginary-real correlations are zero
for f = 0 and f = f_max;
this will make it easier down the road to do matrix interpolation for handling
truncated waveforms.

Although it doesn't matter because of the matrix's symmetry, conceptually (i,j) is (row,column),
and the matrix is stored in column-major format.
(Generally BLAS libraries are optimized for this case.)
*/
class NoiseMatrix
{
 public:

  // Define the index type which is used for indexing matrices at a given frequency.
  // Fastest iteration is over channel; next, over real/imaginary (as unsigned char).
  typedef ProductIndexHandler<RangeIndexHandler<unsigned char>,
                              MapIndexHandler<unsigned char> > NoiseBlockIndexT;

  double& GetCorrByIndex(size_t row, size_t col) {
    assert(row < fNoiseBlockIndex.MaxIndex() and col < fNoiseBlockIndex.MaxIndex());
    return fNoise[row + col*fNoiseBlockIndex.MaxIndex()];
  }

  const double& GetCorrByIndex(size_t row, size_t col) const {
    assert(row < fNoiseBlockIndex.MaxIndex() and col < fNoiseBlockIndex.MaxIndex());
    return fNoise[row + col*fNoiseBlockIndex.MaxIndex()];
  }

  /* Access noise correlations semantically (channel and real/imag). */
  double GetCorrByKey(NoiseBlockIndexT::key_type row,
                      NoiseBlockIndexT::key_type col) const {
    size_t row_index = fNoiseBlockIndex.IndexForKey(row);
    size_t col_index = fNoiseBlockIndex.IndexForKey(col);
    return GetCorrByIndex(row_index, col_index);
  }

 private:
  NoiseMatrix(const NoiseBlockIndexT& index)
    : fNoiseBlockIndex(index),
      fNoise(index.MaxIndex()*index.MaxIndex())
  { }

  NoiseBlockIndexT fNoiseBlockIndex; // Guaranteed equal to the one owned by NoiseCorrelations.
  std::vector<double> fNoise;
  friend class NoiseCorrelations;
};

/*
Class which holds a collection of NoiseMatrices.
This class also holds information about which frequencies are being held, etc.
*/
class NoiseCorrelations
{
 public:

  // Define the index type for frequency information.
  typedef IntervalIndexHandler FrequencyIndexT;

  // Define the index type which is used for indexing matrices at a given frequency.
  // Fastest iteration is over channel; next, over real/imaginary (as unsigned char).
  typedef ProductIndexHandler<RangeIndexHandler<unsigned char>,
                              MapIndexHandler<unsigned char> > NoiseBlockIndexT;

  // To construct an empty noise correlation object, we need to know what channels are desired.
  // Everything else we can set up for you.
  // All correlations are initialized to zero automatically by the resize.
  NoiseCorrelations(const MapIndexHandler<unsigned char>& channelIndex)
    : fFreqIndex(double(1.)/(2048*CLHEP::microsecond), double(0.5)/CLHEP::microsecond, 1024),
      fBlockIndex(RangeIndexHandler<unsigned char>(0, 2), channelIndex)
  {
    fMatrices.resize(fFreqIndex.MaxIndex(), NoiseMatrix(fBlockIndex));
  }

  const FrequencyIndexT& GetFrequencyIndex() const {return fFreqIndex;}

  const NoiseBlockIndexT& GetNoiseBlockIndex() const {return fBlockIndex;}

  NoiseMatrix& GetMatrixForIndex(size_t f)
  {
    assert(f < fFreqIndex.MaxIndex());
    return fMatrices[f];
  }

  const NoiseMatrix& GetMatrixForIndex(size_t f) const
  {
    assert(f < fFreqIndex.MaxIndex());
    return fMatrices[f];
  }

 private:
  FrequencyIndexT fFreqIndex;
  NoiseBlockIndexT fBlockIndex;
  std::vector<NoiseMatrix> fMatrices;
};
#endif
