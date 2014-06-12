#ifndef NoiseMultiplier_hpp
#define NoiseMultiplier_hpp

/*
Class that deals with queueing up vectors to multiply by noise, and performing the
noise multiplication in a multithreaded loop when requested.

At the moment we don't handle interpolating by frequency; but this is where most
of the work to accomplish that would happen.  It would require interpolation
between two noise frequency blocks, which can be done efficiently with the BLAS.
This would be a key component of handling truncated waveforms.

After a call to DoMultiplication, there are three buffers of interest:
fVectorsMultiplied, for the results of multiplication.
fLastInputVectors, the input vectors which were multiplied.
fVectorsToMultiply, the queue for the next multiplication.

These each are really big chunks of memory, and it would be nice to have two instead of three.
This is possible; each event could retain its chunk of memory from fLastInputVectors and,
when done with the last input vector, fill it with the next input vector.
However, the challenge comes from gaps left by events which finish; filling those
gaps efficiently would turn this class into a small memory manager, and although this is
possible (boost has some classes which would help) it is more than I want to take on right now.
Fortunately we aren't too memory-constrained.

fNumVectorsMax is an initial guess for how many columns should be expected; if you set it to a
non-default value, you should do so before calling FillNoiseCorrelations.

FillNoiseCorrelations reads in noise from a file and initializes some data structures.  It should
only be called once.

InsertVector inserts a vector to be multiplied.  A handle is returned, and this should be saved.
You can access the vector in-place before calling DoMultiplication by calling GetInput(size_t handle),
BUT this acquires a shared lock; keep such sections of code short, then release with ReleaseInput().
Calling InsertVector is thread-safe (using internal locking).

DoMultiplication multiplies the input vectors by D^(-1/2)ND^(-1/2).  The results can be retrieved
(and modified in-place) using GetMulResult(size_t handle).  The original input can also be
accessed (and modified) using GetOldInput(size_t handle).
*/

#include "NoiseCorrelations/NoiseCorrelations.hpp"
#include "NoiseCorrelations/NoiseCorrelationsIO.hpp"

#include <boost/thread/mutex.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/lock_guard.hpp>
#include <boost/thread/thread.hpp>
#include <boost/asio/io_service.hpp>
#include <vector>
#include <cassert>

class NoiseMultiplier
{

 public:

  // Typical number of vectors to expect for a multiplication.
  // Set this before calling FillNoiseCorrelations.
  size_t fNumVectorsMax;

  NoiseMultiplier()
  : fNumVectorsMax(500),
    fVectorLength(0),
    fCurrentQueueSize(0),
    fNumVectors(0),
  {}

  // Get the pointer to your input vector from your handle,
  // allowing you to modify it in-place.
  // HOWEVER: if someone tries to insert a vector *and* the insert requires a resize,
  // that would invalidate the pointer.
  // As a result, we acquire a shared_lock; be sure to release it with ReleaseInput
  // when you are done using the pointer.
  double* GetInput(size_t handle) {
    assert(handle % fVectorLength == 0);
    fResizeInQueueMutex.lock_shared();
    assert(handle < fVectorsToMultiply.size());
    return &fVectorsToMultiply[handle];
  }

  void ReleaseInput() {fResizeInQueueMutex.unlock_shared();}

  double* GetOriginalInput(size_t handle) {
    assert(handle < fLastInputVectors.size());
    assert(handle % fVectorLength == 0);
    return &fLastInputVectors[handle];
  }

  double* GetMulResult(size_t handle) {
    assert(handle < fVectorsMultiplied.size());
    assert(handle % fVectorLength == 0);
    return &fVectorsMultiplied[handle];
  }

  // Fill the noise correlations object.
  // Filename is the location to read from.
  // Also allocates initial space based on the hint from fNumVectors (not required but helpful).
  // Only ever call this once.
  void FillNoiseCorrelations(std::string filename, const MapIndexHandler<unsigned char>& channelMap) {
    // Read in noise correlations.
    fNoiseCorr.SetChannelIndex(channelMap);
    NoiseCorrelationsIO::ReadNoiseCorrelations(filename, fNoiseCorr);

    // Calculate how long each vector will be.
    // For all frequencies but the last, the vectors are full-length.
    // For the last frequency, only real coefficients are used.
    assert(fNoiseCorr.GetNoiseBlockIndex().MaxIndex() % 2 == 0);
    fVectorLength = fNoiseCorr.GetFrequencyIndex().MaxIndex() * fNoiseCorr.GetNoiseBlockIndex().MaxIndex();
    fVectorLength -= fNoiseCorr.GetNoiseBlockIndex().MaxIndex()/2;

    // Preallocate space based on the guess of how many we'll need to handle at one time.
    fVectorsToMultiply.resize(fNumVectorsMax*fVectorLength);
    fVectorsMultiplied.resize(fNumVectorsMax*fVectorLength);

    // Precompute D^(1/2) and D^(-1/2) for the full noise matrix.
    // Don't pre-apply it -- we'll need the original matrix for handling truncated waveforms down the road.
    fSqrtDiagNoise.reserve(fVectorLength);
    for(size_t f = 0; f < fNoiseCorr.GetFrequencyIndex().MaxIndex(); f++) {
      size_t BlockSize = (f+1 == fNoiseCorr.GetFrequencyIndex().MaxIndex() ?
                          fNoiseCorr.GetNoiseBlockIndex().MaxIndex() :
                          fNoiseCorr.GetNoiseBlockIndex().MaxIndex()/2);
      for(size_t i = 0; i < BlockSize; i++) {
        fSqrtDiagNoise.push_back(std::sqrt(fNoiseCorr.GetMatrixForIndex(f).GetCorrByIndex(i, i)));
        // In C++11, we should test (in an assert) that fSqrtDiagNoise.back() is normal.
      }
    }
    assert(fSqrtDiagNoise.size() == fVectorLength);
    fInvSqrtDiagNoise.resize(fSqrtDiagNoise.size());
    for(size_t i = 0; i < fSqrtDiagNoise.size(); i++) {
      fInvSqrtDiagNoise[i] = double(1)/fSqrtDiagNoise[i];
    }

    // fNoiseCorr isn't that big memory-wise, so we can afford to make a copy of it
    // and pre-apply fInvSqrtDiagNoise on the left and right.
    // The result, in thesis notation, is D^(-1/2) N D^(-1/2).
    // Note that we *must* retain the original fNoiseCorr for (future) use on truncated waveforms.
    fNoiseCorr_Precon = fNoiseCorr;
    size_t NextDiagIndex = 0;
    for(size_t f = 0; f < fNoiseCorr_Precon.GetFrequencyIndex().MaxIndex(); f++) {
      size_t BlockSize = (f+1 == fNoiseCorr_Precon.GetFrequencyIndex().MaxIndex() ?
                          fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex() :
                          fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex()/2);
      for(size_t i = 0; i < BlockSize; i++) {
        for(size_t j = 0; j < BlockSize; j++) {
          fNoiseCorr_Precon.GetMatrixForIndex(f).GetCorrByIndex(i, j) *= fInvSqrtDiagNoise[NextDiagIndex];
          fNoiseCorr_Precon.GetMatrixForIndex(f).GetCorrByIndex(j, i) *= fInvSqrtDiagNoise[NextDiagIndex];
          NextDiagIndex++;
        }
      }
    }
    assert(NextDiagIndex == fVectorLength);
  }

  // Queue up the data in vec to be multiplied by noise.
  // Calling this function is thread-safe; it locks internally.
  // vec.size() should be a non-negative multiple of the number of rows expected.
  // We return a "handle" (the index where the vector was inserted),
  // which can be used later to retrieve a pointer to the result.
  size_t InsertVector(const std::vector<double>& vec) {
    assert(fVectorLength != 0);
    assert(vec.size() % fVectorLength == 0);
    size_t NumCols = vec.size() / fVectorLength;
    size_t StartIndex;
    assert(NumCols > 0);
    {
      // Reserve space for ourselves in only one thread at a time.
      boost::lock_guard lock(fInsertMutex);
      if(fCurrentQueueSize + vec.size() > fVectorsToMultiply.size()) {
        // Make sure no other threads are using bare pointers to fVectorsToMultiply.
        // This is expensive, but should rarely happen since the size of fVectorsToMultiply only increases.
        boost::lock_guard resizeLock(fResizeInQueueMutex);
        fVectorsToMultiply.resize(fCurrentQueueSize + vec.size());
      }
      StartIndex = fCurrentQueueSize;
      fCurrentQueueSize += vec.size();
      fNumVectors += NumCols;
    }

    // Copy memory (in multithreaded code).
    std::memcpy(reinterpret_cast<void*>(&fVectorsToMultiply[StartIndex]),
                reinterpret_cast<const void*>(&vec[0]),
                sizeof(double)*vec.size());

    return StartIndex;
  }

  // Perform multiplication; start many threads internally.
  // Of course, you should only call this from one thread, and only when you're not doing
  // anything else with the class.
  void DoMultiplication() {
    if(fVectorsMultiplied.size() < fCurrentQueueSize) {
      fVectorsMultiplied.resize(fCurrentQueueSize);
    }

    // Launch threads.
    // Note that newer versions of boost may have a thread pool, which would be easier.
    boost::asio::io_service ioService;
    boost::thread_group thread_group;
    boost::asio::io_service::work* work = new boost::asio::io_service::work(ioService);
    for(size_t i = 0; i < 24; i++) {
      thread_group.create_thread(boost::bind(&boost::asio::io_service::run, &ioService));
    }
    for(size_t f = 0; f < fNoiseCorr_Precon.GetFrequencyIndex().MaxIndex(); f++) {
      ioService.post(boost::bind(&NoiseMultiplier::DoMultiplication_OneFrequency, this, f));
    }
    delete work; // Let ioService know that no more work is coming.
    // FixME: better is to only start 23 threads, since we're the 24th;
    // Then we should (probably) call run from this thread as well.
    // We return when there is no more work for this thread, then join the other threads to verify all are done.
    // It's not a big deal if the join function is not a busy join, but I don't know whether that's the case.
    thread_group.join_all();

    // Prepare to accept more input.
    std::swap(fVectorsToMultiply, fLastInputVectors);
    fCurrentQueueSize = 0;
    fNumVectors = 0;
  }

 private:
  size_t fVectorLength;
  size_t fCurrentQueueSize;
  size_t fNumVectors;
  std::vector<double> fVectorsToMultiply;
  std::vector<double> fVectorsMultiplied;
  std::vector<double> fLastInputVectors;
  NoiseCorrelations fNoiseCorr;
  NoiseCorrelations fNoiseCorr_Precon;
  std::vector<double> fSqrtDiagNoise;
  std::vector<double> fInvSqrtDiagNoise;

  // Mutexes to guarantee thread-safety.
  boost::mutex fInsertMutex;
  boost::shared_mutex fResizeInQueueMutex; // Acquires an exclusive lock before resizing fVectorsToMultiply.

  // Do one multiplication block.
  // This function is thread-safe with no internal locks.
  // f should be a frequency index.
  // Also apply the diagonal preconditioner.
  void DoMultiplication_OneFrequency(size_t f) {
    // We use dgemm -- the noise matrices are symmetric, but often dsymm is slower than dgemm.
    // But, we should test this to verify that dsymm isn't faster than dgemm.
    assert(fVectorsToMultiply.size() == fVectorsMultiplied.size());
    assert(fVectorsToMultiply.size() == fNumVectors*fVectorLength);
    assert(fInvSqrtDiagNoise.size() == fVectorLength);

    // Find the address of the first relevant entry in fVectorsToMultiply and fVectorsMultiplied.
    size_t VectorIndex = f*fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex();
    double* InVectorLoc = &fVectorsToMultiply[VectorIndex];
    double* OutVectorLoc = &fVectorsMultiplied[VectorIndex];

    // Find the height of this block.
    size_t BlockSize = (f+1 == fNoiseCorr_Precon.GetFrequencyIndex().MaxIndex() ?
                        fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex() :
                        fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex()/2);
    assert(VectorIndex + BlockSize < fVectorLength);

    // Pointer to the first relevant entry of the noise matrix block.
    const double* NoiseMat = &fNoiseCorr_Precon.GetMatrixByIndex(f).GetCorrByIndex(0, 0);

    // Span of noise block.
    size_t NoiseSpan = fNoiseCorr_Precon.GetNoiseBlockIndex().MaxIndex();

    // Check that conversions from size_t to MKL_INT are valid.
    assert(BlockSize     <= std::numeric_limits<MKL_INT>::max() and
           fNumVectors   <= std::numeric_limits<MKL_INT>::max() and
           NoiseSpan     <= std::numeric_limits<MKL_INT>::max() and
           fVectorLength <= std::numeric_limits<MKL_INT>::max());

    // Do the multiplication.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans
                BlockSize, fNumVectors, BlockSize,
                double(1), NoiseMat, NoiseSpan,
                InVectorLoc, fVectorLength,
                double(0), OutVectorLoc, fVectorLength);
  }
};
#endif
