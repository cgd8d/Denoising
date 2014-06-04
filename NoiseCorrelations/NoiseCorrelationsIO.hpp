#ifndef NoiseCorrelationsIO_hpp
#define NoiseCorrelationsIO_hpp

/*
Read a noise correlations object from memory, pack it according to
available symmetries (which are not exploited in memory due to speed concerns),
and write to disk.


NOTE: comments below are just me planning things out.

I am choosing to use HDF5.  I was going to use boost.serialization, but turns out that doesn't provide a portable binary format, so it's not really addressing my problem.  ROOT is definitely not the right thing for storing big raw-array datasets like this.

I will keep this very simple: in the "/" group (root group) I store datasets called noise_correlations_<f>, where f is the frequency index of the dataset.  <f> should range from 0 to N-1.  Additionally, and also in the root directory, I will store a vector of channel keys to indicate the channel map (called "channel_map").  Finally, the root group will have attributes: f_min (floating-point), f_max (floating-point), and num_f (unsigned integral), indicating how many frequency datasets there are and their semantic meaning.  (The frequency corresponding to dataset i is f_min + (f_max-f_min)*i/num_f, so f_max is non-inclusive, just as in the index handler.)
*/

#include "NoiseCorrelations/NoiseCorrelations.hpp"
#include "Utilities/IndexHandler.hpp"
#include "hdf5.h"
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

/*
When an error occurs in HDF5, print a message, then die;
I don't want to have to handle failures.
This allows me to avoid explicitly checking return values throughout.
*/
extern "C"
herr_t custom_HDF5_error_callback(hid_t estack_id, void* /* unused */)
{
  H5Eprint2(estack_id, NULL);
  exit(1);
  return 1;
}

namespace NoiseCorrelationIO
{

/*
Initialize HDF5 at program initialization to use custom_HDF5_error_callback for errors.
*/
struct custom_HDF5_init {
  custom_HDF5_init() {
    herr_t ret = H5Eset_auto2(H5E_DEFAULT, &custom_HDF5_error_callback, NULL);
    if(ret < 0) {
      std::cerr << "Failed to set HDF5 error callback." << std::endl;
      exit(1);
    }
  }
} g_custom_HDF5_init;


/*
Write an existing noise correlation object to disk.
If filename is already a file, fail.
*/
void WriteNoiseCorrelations(std::string filename, const NoiseCorrelations& noise)
{
  hid_t fileID = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

  // Write the channel index as an ordered list of included channels.
  {
    const MapIndexHandler<unsigned char>& chanMap = noise.GetNoiseBlockIndex().MinorIndex();
    std::vector<unsigned char> chanVec; // Make a buffer to write.
    for(size_t i = 0; i < chanMap.MaxIndex(); i++) chanVec.push_back(chanMap.KeyAtIndex(i));
    const hsize_t NumChannels = chanMap.MaxIndex();
    hid_t vectorID = H5Screate_simple(1, &NumChannels, NULL);
    hid_t attID = H5Acreate2(fileID, "channel_list", H5T_STD_U8LE, vectorID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(&chanVec[0]));
    H5Aclose(attID);
    H5Sclose(vectorID);
  }

  // Write the actual noise information.
  // We choose to write one dataset per frequency, since those are stored in memory as separate arrays.
  // However, we need to manually pick out only the entries which contain non-redundant information.
  std::vector<double> PackedArray; // Reuse rather than re-allocating each time.
  for(size_t f = 0; f < 1024; f++) {
    const NoiseMatrix& mat = noise.GetMatrixForIndex(f);
    const NoiseCorrelations::NoiseBlockIndexT& NoiseBlockIndex = noise.GetNoiseBlockIndex();
    assert(NoiseBlockIndex.MaxIndex() % 2 == 0 and NoiseBlockIndex.MaxIndex() > 0);

    // Create the name for this dataset.
    std::ostringstream strstream;
    strstream << "/noise_corr_" << std::setfill('0') << std::setw(4) << f;
    std::string dataset_name = strstream.str();

    // Allocate space in the temporary vector.
    hsize_t ExpectedSize;
    if(f < 1023) ExpectedSize = NoiseBlockIndex.MaxIndex()*NoiseBlockIndex.MaxIndex()/4;
    else ExpectedSize = NoiseBlockIndex.MaxIndex()*(NoiseBlockIndex.MaxIndex()+2)/8;
    PackedArray.resize(0);
    PackedArray.reserve(NoiseBlockIndex.MaxIndex()*NoiseBlockIndex.MaxIndex()/4);

    // Fill PackedArray. Take into account all appropriate symmetries.
    for(size_t i = 0; i < NoiseBlockIndex.MaxIndex(); i++) {
      for(size_t j = i; j < NoiseBlockIndex.MaxIndex(); j++) {
        if(f < 1023) { // Symmetries for most frequency indices.
          // real-real = imag-imag, so skip imag-imag.
          if(NoiseBlockIndex.MajorIndexForIndex(i) == 1 and
             NoiseBlockIndex.MajorIndexForIndex(j) == 1) continue;
          // real-imag block is antisymmetric.
          if(NoiseBlockIndex.MajorIndexForIndex(i) == 0 and
             NoiseBlockIndex.MajorIndexForIndex(j) == 1 and
             NoiseBlockIndex.MinorIndexForIndex(i) >= NoiseBlockIndex.MinorIndexForIndex(j)) continue;
        }
        else { // Symmetries for the last frequency bin.
          // All imaginary FFT components are zero.
          if(NoiseBlockIndex.MajorIndexForIndex(i) == 1 or
             NoiseBlockIndex.MajorIndexForIndex(j) == 1) continue;
        }

        // If we reach this point, this is a non-redundant noise component.
        PackedArray.push_back(mat.GetCorrByIndex(i, j));
      }
    }
    assert(PackedArray.size() == ExpectedSize);

    // Write the array to file.
    hid_t vectorID = H5Screate_simple(1, &ExpectedSize, NULL);
    hid_t datasetID = H5Dcreate2(fileID, dataset_name.c_str(), H5T_IEEE_F64LE, vectorID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             reinterpret_cast<void*>(&PackedArray[0]));
    H5Dclose(datasetID);
    H5Sclose(vectorID);
  }

  H5Fclose(fileID);
}





} // namespace NoiseCorrelationIO
#endif
