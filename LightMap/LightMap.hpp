#ifndef LightMap_hpp
#define LightMap_hpp

/*
Classes for storing information about the lightmap (position-dependent and time-dependent parts).

There are a couple of helper classes in here, so we put them all in a namespace LightMap.
The user should be aware of:

PositionFunc: Holds R(\vec{x}, apdNo) in one massive array.
GainSnapshot: Holds S(t, apdNo) for one particular range of runs.
GainFunc: Holds S(t, apdNo) for all times as a non-overlapping set of GainSnapshots.

The gain function is assumed to be constant on run intervals [first_run, last_run]; note
that the old version of this code (used for the 2014 0nu paper) assumed gain could change
on a run-by-run basis, but this led to a very jump gain function which probably was
dominated by statistics.
*/

#include "Utilities/IndexHandler.hpp"
#include <vector>
#include <algorithm>
#include <cassert>

namespace LightMap {

const unsigned char MAX_APDS = 71; // The maximum number of APDs we have ever had or ever will have.

// Map indexing the actual APDs present.
typedef MapIndexHandler<unsigned char> APDIndexT;

/*
The lightmap functions are functions of APD number and some other parameter: R(\vec{x}, i) and S(t, i).
Generally we get an event at \vec{x}, t, and want to retrieve R and S for all i at that point.
This is most efficient if function values for different i are stored in adjacent memory locations.
This structure the values of a function for all APDs; it ensures they are adjacent by storing
them in a statically-allocated c-array.  (Note that a std::vector would be *wrong* for this,
since it only contains a pointer; the data is allocated from elsewhere.)
Since I don't know the number of APDs we're using, I have to make the array as big as the maximum
possible number of APDs in use, but generally that will be close to the true number.
In C++11, std::array would serve this purpose.
Note: don't store the true number of APDs here, the overhead is unreasonable.  Store it in
the container classes via an index.
*/
struct FuncVsAPD
{
  double fVals[MAX_APDS];
  double& operator[](size_t i) {return fVals[i];}
  const double& operator[](size_t i) const {return fVals[i];}
};

/*
Container for the function which translates from position to lightmap yield.
*/
class PositionFunc
{
 public:
  // Map three-dimensional points to indices.
  typedef ProductIndexHandler<ProductIndexHandler<IntervalIndexHandler, IntervalIndexHandler>,
                              IntervalIndexHandler> PosIndexT;

  // When you construct the position map, you should already know what binning you're
  // using and which APDs are present.
  PositionFunc()
  : fPosIndex(XYPosIndexT(IntervalIndexHandler(0, 1, 1), IntervalIndexHandler(0, 1, 1)),
              IntervalIndexHandler(0, 1, 1)),
    fIsInitialized(false)
  { }

  // Position binning to use, in units of mm.
  // Shouldn't change this once the map is filled.
  void SetBinning(double xmin, double xstep, size_t nx,
                  double ymin, double ystep, size_t ny,
                  double zmin, double zstep, size_t nz) {
    fPosIndex = PosIndexT(XYPosIndexT(IntervalIndexHandler(xmin, xstep, nx),
                                      IntervalIndexHandler(ymin, ystep, ny)),
                          IntervalIndexHandler(zmin, zstep, nz));
    fData.resize(fPosIndex.MaxIndex());
    fIsInitialized = true;
  }

  void SetAPDIndex(const APDIndexT& index) {
    assert(index.MaxIndex() <= MAX_APDS);
    fAPDIndex = index;
  }

  double& GetValAt(size_t pos_index, size_t apd_index) {
    assert(pos_index < fPosIndex.MaxIndex());
    assert(apd_index < fAPDIndex.MaxIndex());
    assert(fIsInitialized);
    return fData[pos_index][apd_index];
  }

  const double& GetValAt(size_t pos_index, size_t apd_index) const {
    assert(pos_index < fPosIndex.MaxIndex());
    assert(apd_index < fAPDIndex.MaxIndex());
    assert(fIsInitialized);
    return fData[pos_index][apd_index];
  }

  FuncVsAPD& GetAllValsAt(size_t pos_index) {
    assert(pos_index < fPosIndex.MaxIndex());
    assert(fIsInitialized);
    return fData[pos_index];
  }

  const FuncVsAPD& GetAllValsAt(size_t pos_index) const {
    assert(pos_index < fPosIndex.MaxIndex());
    assert(fIsInitialized);
    return fData[pos_index];
  }

  const PosIndexT& PosIndex() const {
    assert(fIsInitialized);
    return fPosIndex;
  }

  const APDIndexT& APDIndex() const {return fAPDIndex;}

 private:
  typedef ProductIndexHandler<IntervalIndexHandler, IntervalIndexHandler> XYPosIndexT;
  PosIndexT fPosIndex;
  APDIndexT fAPDIndex;
  bool fIsInitialized;
  std::vector<FuncVsAPD> fData;
};

/*
Store the gain function for a particular range of time.
This just consists of a FuncVsAPD plus the metadata to identify the time of validity etc.
Since there aren't too many gain snapshots in any one program, it doesn't need to be compact.
The time of validity is from fFirstRun to fLastRun, inclusive.
*/
class GainSnapshot
{
 public:

  GainSnapshot(APDIndexT index, int first_run, int last_run)
  : fAPDIndex(index),
    fFirstRun(first_run),
    fLastRun(last_run)
  {
    assert(fAPDIndex.MaxIndex() <= MAX_APDS and first_run <= last_run);
  }

  // Implement a partial ordering of gain snapshots.
  bool operator<(const GainSnapshot& other) const {return fLastRun < other.fFirstRun;}
  bool operator>(const GainSnapshot& other) const {return fFirstRun > other.fLastRun;}

  // Also implement a partial ordering with respect to run numbers.
  bool operator<(const int& run) const {return fLastRun < run;}
  bool operator>(const int& run) const {return fFirstRun > run;}

  double& GetValAt(size_t apd_index) {
    assert(apd_index < fAPDIndex.MaxIndex());
    return fData[apd_index];
  }

  const double& GetValAt(size_t apd_index) const {
    assert(apd_index < fAPDIndex.MaxIndex());
    return fData[apd_index];
  }

  const APDIndexT& APDIndex() const {return fAPDIndex;}

  int FirstRun() const {return fFirstRun;}

  int LastRun() const {return fLastRun;}

 private:
  APDIndexT fAPDIndex;
  FuncVsAPD fData;
  int fFirstRun;
  int fLastRun;
};

/*
In cases where we need to store gain functions covering all history, store them here.
*/
class GainFunc
{
 public:
  GainFunc() { }

  // Snapshots must not overlap, but this is checked when the function is used.
  // All snapshots should have APD indices with identical sets of keys,
  // though the ordering may differ.
  void InsertGain(const GainSnapshot& snapshot) {
    if(NumSnapshots() > 0) {
      // In assertions, verify that the APD index is consistent with what we already have.
      assert(snapshot.APDIndex().MaxIndex() == fSnapshots[0].APDIndex().MaxIndex());
      for(size_t i = 0; i < fSnapshots[0].APDIndex().MaxIndex(); i++) {
        assert(snapshot.APDIndex().HasKey(fSnapshots[0].APDIndex().KeyForIndex(i)));
      }
    }
    fSnapshots.push_back(snapshot);
    Sort(); // Inefficient to call this repeatedly, but generally the effect is negligible.
  }

  // Convenience function to make a temporary snapshot and insert it.
  void InsertGain(const APDIndexT& index, int first_run, int last_run) {
    GainSnapshot gain(index, first_run, last_run);
    InsertGain(gain);
  }

  size_t NumSnapshots() const {return fSnapshots.size();}

  GainSnapshot& GainAtIndex(size_t index) {
    assert(index < fSnapshots.size());
    return fSnapshots[index];
  }

  const GainSnapshot& GainAtIndex(size_t index) const {
    assert(index < fSnapshots.size());
    return fSnapshots[index];
  }

  GainSnapshot& GainForRun(int run) {
    std::vector<GainSnapshot>::iterator it =
      std::lower_bound(fSnapshots.begin(), fSnapshots.end(), run);
    assert(it != fSnapshots.end());
    assert(not (*it > run)); // Verify that it really is in the range.
    return *it;
  }

  const GainSnapshot& GainForRun(int run) const {
    std::vector<GainSnapshot>::const_iterator it =
      std::lower_bound(fSnapshots.begin(), fSnapshots.end(), run);
    assert(it != fSnapshots.end());
    assert(not (*it > run)); // Verify that it really is in the range.
    return *it;
  }

 private:
  std::vector<GainSnapshot> fSnapshots;

  void Sort() {
    std::sort(fSnapshots.begin(), fSnapshots.end());
    for(size_t i = 1; i < fSnapshots.size(); i++) {
      assert(fSnapshots[i-1] < fSnapshots[i]); // Ensure there are no overlaps.
    }
  }
};

} // namespace LightMap
#endif
