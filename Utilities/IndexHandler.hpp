/*
Manage the association with indices (a consecutive numbering scheme 0,1,...,max-1)
to a key.  The index indicates memory layout, and the key describes
semantical meaning.

This is currently managed with three template classes:
MapIndexHandler:	index a discrete set of keys.
IntervalIndexHandler:	index evenly spaced values from a floating-point interval.
RangeIndexHandler:	index evenly spaced values from a range of integers.
ProductIndexHandler:	index the cartesian product of two sets of keys.
SubIndexHandler:	index which uses only a subrange of another index.

All of these classes expose a common interface:
key_type:				The type of the keys (a typedef).
key_type KeyForIndex(size_t) const:	Convert index to key.
size_t MaxIndex() const:		Index ranges from 0 to MaxIndex()-1, inclusive.

Additional functions specific to the type of handler is also provided.

In principle this common interface could be exposed by a virtual base class,
but these are all meant to be *very* lightweight classes to simplify the management
of semantics for very cheap operations; adding a virtual interface would mean
lookup through a vtable, which would be overkill.
*/

#ifndef IndexHandler_hpp
#define IndexHandler_hpp

#include <cassert>
#include <limits>
#include <vector>
#include <algorithm>

/*
Manage conversion from some arbitrary set of integral values (eg. channel number)
to indices.

Keys should be inserted one by one; they will be assigned indices in order.

KeyForIndex is currently inefficient; it could be faster if we maintained a secondary
structure which maps back.  But it complicates the code and I may not need it, so omit for now.
*/
template<typename KeyT>
class MapIndexHandler
{
 public:

  typedef KeyT key_type;

  MapIndexHandler() : fMaxIndex(0)
  {
    assert(std::numeric_limits<key_type>::is_integer); // FixME: make this a compile-time test.
  }

  // Add a new key.  Duplicate keys should not be added.
  void InsertKey(const key_type& k) {
    assert(std::find(fKeys.begin(), fKeys.end(), k) == fKeys.end());
    fKeys.push_back(k);
    fMaxIndex++;
  }

  key_type KeyAtIndex(size_t index) const {
    assert(index < fKeys.size());
    return fKeys[index];
  }

  bool HasKey(key_type key) const {
    return (std::find(fKeys.begin(), fKeys.end(), key) != fKeys.end());
  }

  // Not efficient, but simplest from a coding perspective -- be sure to profile.
  // (If it is a bottleneck, consider a secondary structure or boost::MultiIndex.)
  // Do not pass in a key which is not part of the index.
  size_t IndexForKey(key_type key) const {
    for(size_t i = 0; i < fMaxIndex; i++) {
      if(fKeys[i] == key) return i;
    }

    assert(false and "key is not in index"); // Should never reach this point.
    return 0; // Just to make the compiler happy.
  }

  size_t MaxIndex() const {return fMaxIndex;}

 private:
  size_t fMaxIndex;
  std::vector<KeyT> fKeys;
};

/*
For an interval [a,b] and N indices, translate index:
i -> a + (b-a)*i/(N-1), for i = 0 .. N-1.
This means we divide the interval into N even pieces and index including a and b.
Mainly this is useful for indexing frequency and keeping track of the meanings of frequency.
*/
class IntervalIndexHandler
{
 public:
  typedef double key_type;

  IntervalIndexHandler(double start, double end, size_t NumIndices)
    : fMaxIndex(NumIndices),
      fStart(start),
      fStepSize((end-start)/NumIndices)
  {
    assert(fStepSize > 0 and fMaxIndex > 0);
  }

  size_t MaxIndex() const {return fMaxIndex;}

  key_type KeyForIndex(size_t index) const {
    assert(index < fMaxIndex);
    return fStart + fStepSize*index;
  }

  // Return the largest index such that KeyForIndex(index) <= key.
  // If key < start, return size_t(-1).
  size_t IndexForKey(key_type key) const {
    key -= fStart;
    key /= fStepSize;
    if(key < 0) return size_t(-1); // Defer to here to avoid floating-point issues.
    size_t entry = size_t(key);
    if(entry >= fMaxIndex) return fMaxIndex-1;
    return entry;
  }

  double Start() const {return fStart;}
  double End() const {return fStart + fMaxIndex*fStepSize;}
  double StepSize() const {return fStepSize;}

 private:
  size_t fMaxIndex; // Should this be made a template argument?
  double fStart;
  double fStepSize;
};

/*
Translate a range of integral values k1 ... kN to indices 0 ... N-1.
The keys must be consecutive -- we don't allow strides.  (This could be relaxed.)
Specify the range with Start and NumKeys, where End is not inclusive.
FixME: does allowing step sizes other than one hurt speed?  And is it useful?
*/
template<typename KeyT>
class RangeIndexHandler
{
 public:
  typedef KeyT key_type;

  RangeIndexHandler(KeyT start, size_t NumKeys)
    : fMaxIndex(NumKeys),
      fStart(start)
  {
    assert(NumKeys > 0);

    // Verify that we're not exceeding the range of the key_type.
    // This is surprisingly tricky -- beware the cases where key_type is signed or bool.
    assert(size_t(std::numeric_limits<key_type>::max()) - start + 1 >= NumKeys);
  }

  size_t MaxIndex() const {return fMaxIndex;}

  key_type KeyForIndex(size_t index) const {
    assert(index < fMaxIndex);
    return fStart + key_type(index);
  }

  size_t IndexForKey(key_type key) const {
    assert(key >= fStart);
    key -= fStart;
    size_t ret = size_t(key);
    assert(ret < fMaxIndex);
    return ret;
  }

 private:
  size_t fMaxIndex;
  key_type fStart;
};

/*
Often it may be the case that we wish to index more than one key.
This class takes in two existing index handlers and uses
the first as a major index and the second as a minor index,
so the keys for the combination is the cartesian product of
keys for each individual index.
*/
template<class MajorIndexT, class MinorIndexT>
class ProductIndexHandler
{
 public:
  typedef std::pair<typename MajorIndexT::key_type,
                    typename MinorIndexT::key_type> key_type;

  ProductIndexHandler(const MajorIndexT& index1, const MinorIndexT& index2) :
    fMaxIndex(index1.MaxIndex()*index2.MaxIndex()),
    fMajorIndex(index1),
    fMinorIndex(index2)
  { }

  size_t MaxIndex() const {return fMaxIndex;}

  key_type KeyForIndex(size_t index) const {
    assert(index < fMaxIndex);
    size_t i1 = MajorIndexForIndex(index);
    size_t i2 = MinorIndexForIndex(index);
    return std::make_pair(fMajorIndex.KeyForIndex(i1), fMinorIndex.KeyForIndex(i2));
  }

  size_t IndexForKey(key_type key) const {
    size_t k1 = fMajorIndex.IndexForKey(key.first);
    size_t k2 = fMinorIndex.IndexForKey(key.second);
    return k1*fMinorIndex.MaxIndex() + k2;
  }

  const MajorIndexT& MajorIndex() const {return fMajorIndex;}
  size_t MajorIndexForIndex(size_t index) const {
    return index/fMinorIndex.MaxIndex();
  }

  const MinorIndexT& MinorIndex() const {return fMinorIndex;}
  size_t MinorIndexForIndex(size_t index) const {
    return index%fMinorIndex.MaxIndex();
  }

 private:
  size_t fMaxIndex;
  MajorIndexT fMajorIndex;
  MinorIndexT fMinorIndex;
};

/*
Use a sub-range [start, end) of a given index.
*/
template<class IndexT>
class SubIndexHandler
{
 public:
  typedef typename IndexT::key_type key_type;

  SubIndexHandler(const IndexT& index, size_t start, size_t end) :
    fMaxIndex(end-start),
    fStart(start),
    fBaseIndex(index)
  {
    assert(start < end);
  }

  size_t MaxIndex() const {return fMaxIndex;}

  key_type KeyForIndex(size_t index) const {
    assert(index < fMaxIndex);
    return fBaseIndex.KeyForIndex(index + fStart);
  }

  size_t IndexForKey(key_type key) const {
    size_t index = fBaseIndex.IndexForKey(key);
    assert(fStart <= index and index < fStart + fMaxIndex);
    return index - fStart;
  }

  const IndexT& BaseIndex() const {return fBaseIndex;}

 private:
  size_t fMaxIndex;
  size_t fStart;
  IndexT fBaseIndex;
};
#endif
