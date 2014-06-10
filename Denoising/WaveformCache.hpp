#ifndef WaveformCache_hpp
#define WaveformCache_hpp

/*
Class for storing waveforms which will be required later.  Waveforms need to be used at
two separate points in processing -- once for reconstruction, and once for applying
denoising.  We cache the waveforms rather than re-reading them from disk.

We assume that waveforms will only need to be retrieved once; therefore, the semantics are:
1) void SetChannelIndex(const MapIndexHandler<unsigned char>& index);
2) void CacheWaveforms(const EXOEventData* event);
3) WaveformCache::WaveformHolderT* GetWaveforms(int runNo, int eventNo);
You are responsible for deleting the object returned by GetWaveforms.

Currently we store the waveforms as unsigned shorts, and only store the waveforms
which we'll need; this seems to make them compact enough that we don't exhaust memory.
However, if that changed (eg. because it took more iterations to denoise), there are
a few forms of compression we could try:

1) Simplest would be to squeeze two samples into 3 bytes rather than 4; since ADC counts are
taken as 12-bit values, this works.  This isn't very much compression, though, only 1.3x

2) zlib is available on all computing systems we'd care about, and JJ has estimated it provides
a factor of 4.3x compression compared to unsigned shorts.  Very available, very well-tested,
and can make use of multi-threading.

3) JJ has the compression algorithm which is used in online code; alternatively, he has suggested
that it would be possible to grab the compressed data directly (since it is initially read in
compressed form), which would save the processing associated with compression as well.  This gives
the best compression factor (6.6x, compared to unsigned shorts) and is the fastest, but would
require developer time.

4) In principle we could use EXOWaveform::Compress and Decompress; however, I have some concerns about
the portability of that code and how carefully tested it is, and I would suggest using something more
proven.  (Why?  Well, for example I noticed bitwise operations being performed on signed values,
which is not a portable operation.  Probably it works fine, but in general I think we do too much
in-house which needn't be done in-house at all.)  If you want to go this route, though, go for it;
I'm not in the collaboration anymore.
*/

#include "Utilities/IndexHandler.hpp"
#include "EXOUtilities/EXOEventData.hh"
#include "EXOUtilities/EXOWaveformData.hh"
#include "EXOUtilities/EXOWaveform.hh"
#include <map>
#include <vector>
#include <utility>
#include <cassert>

class WaveformCache
{
 public:

  // Hold cached waveforms for a single event.
  typedef std::vector<std::vector<unsigned short> > WaveformHolderT;

  void SetChannelIndex(const MapIndexHandler<unsigned char>& index) {fChannelIndex = index;}
  const MapIndexHandler<unsigned char>& ChannelIndex() const {return fChannelIndex;}

  void CacheWaveforms(const EXOEventData* event) {
    assert(event);
    int runNo = event->fRunNumber;
    int eventNo = event->fEventNumber;

    WaveformHolderT* cachedWFs = new WaveformHolderT;
    cachedWFs.resize(fChannelIndex.MaxIndex());
    for(size_t i = 0; i < fChannelIndex.MaxIndex(); i++) {
      const EXOWaveform* wf = event->GetWaveformData()->GetWaveformWithChannel(fChannelIndex.KeyForIndex(i));
      assert(wf->GetLength() == event->GetWaveformData()->fNumSamples);
      (*cachedWFs)[i].resize(event->GetWaveformData()->fNumSamples);
      for(size_t t = 0; t < wf->GetLength(); t++) {
        assert(wf[t] >= 0 and wf[t] < 4096);
        (*cachedWFs)[i][t] = (unsigned short)wf[t];
      }
    }

    std::pair<int, int> key(runNo, eventNo);
    assert(fCachedWaveforms.find(key) == fCachedWaveforms.end());
    fCachedWaveforms[key] = cachedWFs;
  }

  // Once you call this function, *you* own the returned pointer.
  // Retrieving a run/event number which have not been cached is a fatal error.
  WaveformHolderT* GetWaveforms(int runNo, int eventNo) {
    CacheContainerT::iterator it = fCachedWaveforms.find(std::make_pair(runNo, eventNo));
    assert(it != fCachedWaveforms.end());
    WaveformHolderT* ret = it->second;
    fCachedWaveforms.erase(it);
    return ret;
  }

  ~WaveformCache() {
    for(CacheContainerT::iterator it = fCachedWaveforms.begin();
        it != fCachedWaveforms.end();
        it++) {
      delete it->second;
    }
  }

 private:
  typedef std::map<std::pair<int, int>, WaveformHolderT*> CacheContainerT;
  MapIndexHandler<unsigned char> fChannelIndex;
  CacheContainerT fCachedWaveforms;
};
#endif
