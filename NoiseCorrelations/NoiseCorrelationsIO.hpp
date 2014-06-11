#ifndef NoiseCorrelationsIO_hpp
#define NoiseCorrelationsIO_hpp

/*
Handle reading and writing NoiseCorrelations objects.  These functions take care of
packing the noise information to be more compact on disk.

Users will generally only call the functions:

namespace NoiseCorrelationsIO {
  void WriteNoiseCorrelations(std::string filename, const NoiseCorrelations& noise);
  void ReadNoiseCorrelations(std::string filename, NoiseCorrelations& noise);
}

In each case, the in-memory channel map is read from "noise".  WriteNoiseCorrelations
will fail if the file already exists.  WriteNoiseCorrelations expects the noise
matrices to obey all of the expected symmetries, and as a result it does not examine
entries which are expected to be redundant; if the symmetries are not obeyed, then
the matrix on disk will not reflect the version in memory.

I currently always assume that the noise is computed from 2048-sample waveforms,
and that we save all except the 0-frequency component.  In principle it would be nice
to generalize that, but in practice it is helpful to guarantee that for full waveforms
the frequencies in NoiseCorrelations exactly correspond to the frequencies in data;
otherwise we end up doing more noise interpolation than is necessary.
But this will bite us if we later experiment with different-length waveforms.

We use HDF5 because of its strong guarantees of binary-format floating-point portability
and good performance.  (Note: boost.serialization does not guarantee floating-point
portability in binary format, as of boost 1.55.  ROOT guarantees portability,
but its performance isn't optimized for large arrays.)

The root directory contains:
A version number (currently 1) in case the format is later modified.
A channel map to indicate which channels have noise correlations saved.
One array per frequency, representing the noise matrix, packed to account for symmetries.

All errors in HDF5 cause termination of the program, with an accompanying error message.
*/

#include "NoiseCorrelations/NoiseCorrelations.hpp"
#include "Utilities/IndexHandler.hpp"
#include "Utilities/HDF5Helper.hpp"
#include "hdf5.h"
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cassert>

namespace NoiseCorrelationIO
{

size_t ExpectedPackedSize(size_t f, size_t UnpackedSize) {
  assert(UnpackedSize % 2 == 0 and UnpackedSize > 0);
  if(f < 1023) return UnpackedSize*UnpackedSize/4;
  else return UnpackedSize*(UnpackedSize+2)/8;
}

bool IncludeEntryInPackedArray(size_t f, size_t row, size_t col,
                               const NoiseCorrelations::NoiseBlockIndexT& index) {
  assert(row < index.MaxIndex() and col < index.MaxIndex() and index.MaxIndex() % 2 == 0);
  if(row > col) return false;

  if(f < 1023) { // Symmetries for most frequency indices.
    // real-real = imag-imag, so skip imag-imag.
    if(index.MajorIndexForIndex(row) == 1 and
       index.MajorIndexForIndex(col) == 1) return false;
    // real-imag block is antisymmetric.
    if(index.MajorIndexForIndex(row) == 0 and
       index.MajorIndexForIndex(col) == 1 and
       index.MinorIndexForIndex(row) >= index.MinorIndexForIndex(col)) return false;
  }

  else { // Symmetries for the last frequency bin.
    // All imaginary FFT components are zero.
    if(index.MajorIndexForIndex(row) == 1 or
       index.MajorIndexForIndex(col) == 1) return false;
  }

  return true;
}

// Return the index in the packed array where we should get this entry.
// ret.first is the index; ret.second is true if we'll need a sign flip.
// ret.first == size_t(-1) if this entry is identically zero and not stored.
std::pair<size_t, bool> PackedArrayIndexFor(size_t f, size_t row, size_t col,
                                            const NoiseCorrelations::NoiseBlockIndexT& index) {
  assert(row < index.MaxIndex() and col < index.MaxIndex() and index.MaxIndex() % 2 == 0);
  if(row > col) std::swap(row, col);
  bool flip_sign = false;

  if(f < 1023) {
    if(index.MajorIndexForIndex(row) == 1 and
       index.MajorIndexForIndex(col) == 1) {
      // Given imag-imag, retrieve it from real-real instead.
      row = row % index.MinorIndex().MaxIndex();
      col = col % index.MinorIndex().MaxIndex();
    }
    else if(index.MajorIndexForIndex(row) == 0 and
            index.MajorIndexForIndex(col) == 1) {
      // We're in the real-imag block, which is antisymmetric.
      if(index.MinorIndexForIndex(row) == index.MinorIndexForIndex(col)) {
        return std::make_pair(size_t(-1), flip_sign);
      }
      if(index.MinorIndexForIndex(row) > index.MinorIndexForIndex(col)) {
        size_t tmp = row;
        row = col % index.MinorIndex().MaxIndex();
        col = index.MinorIndex().MaxIndex() + tmp % index.MinorIndex().MaxIndex();
        flip_sign = true;
      }
    }

    // Compute index.
    size_t packed_index = row*(2*index.MinorIndex().MaxIndex() - row);
    packed_index += col - row;
    if(index.MajorIndexForIndex(col) == 1) packed_index -= row + 1;
    return std::make_pair(packed_index, flip_sign);
  }

  else { // f == 1023
    if(index.MajorIndexForIndex(row) == 1 or
       index.MajorIndexForIndex(col) == 1) {
         // All imaginary components are identically zero.
         return std::make_pair(size_t(-1), flip_sign);
    }

    // That's the only identity; compute index.
    size_t packed_index = row*index.MinorIndex().MaxIndex() - row*(row-1)/2;
    packed_index += col - row;
    return std::make_pair(packed_index, flip_sign);
  }
}

/*
Write an existing noise correlation object to disk.
If filename is already a file, fail.
*/
void WriteNoiseCorrelations(std::string filename, const NoiseCorrelations& noise)
{
  hid_t fileID = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

  // Write a version number of 1.
  {
    unsigned char version = 1;
    hid_t scalarID = H5Screate(H5S_SCALAR);
    hid_t attID = H5Acreate2(fileID, "version", H5T_STD_U8LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(&version));
    H5Aclose(attID);
    H5Sclose(scalarID);
  }

  // Write the channel index as an ordered list of included channels.
  HDF5Helper::WriteMapAsAttribute(noise.GetNoiseBlockIndex().MinorIndex(),
                                  fileID,
                                  "channel_list");

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
    hsize_t ExpectedSize = ExpectedPackedSize(f, NoiseBlockIndex.MaxIndex());
    PackedArray.resize(0);
    PackedArray.reserve(NoiseBlockIndex.MaxIndex()*NoiseBlockIndex.MaxIndex()/4);

    // Fill PackedArray. Take into account all appropriate symmetries.
    for(size_t i = 0; i < NoiseBlockIndex.MaxIndex(); i++) {
      for(size_t j = i; j < NoiseBlockIndex.MaxIndex(); j++) {
        if(not IncludeEntryInPackedArray(f, i, j, NoiseBlockIndex)) continue;
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

  assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
  H5Fclose(fileID);
}

/*
Read noise correlations from a file into a noise object.
Only read noise related to channels which are in the noise object's map,
so be sure to call SetChannelIndex before reading correlations into it.
*/
void ReadNoiseCorrelations(std::string filename, NoiseCorrelations& noise)
{
  hid_t fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Verify this is version 1.
  {
    unsigned char version;
    hid_t attID = H5Aopen(fileID, "version", H5P_DEFAULT);
    H5Aread(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(version));
    H5Aclose(attID);
    assert(version == 1);
  }

  // I have one index for the in-memory noise matrices; build another for the file.
  MapIndexHandler<unsigned char> FileChannelMap = HDF5Helper::ReadMapFromAttribute(fileID, "channel_list");
  for(size_t i = 0; i < noise.GetNoiseBlockIndex().MaxIndex(); i++) {
    // Verify that the file has correlation information for every channel we're requesting.
    assert(FileChannelMap.HasKey(noise.GetNoiseBlockIndex().KeyForIndex(i).second));
  }
  NoiseCorrelations::NoiseBlockIndexT FileBlockIndex(RangeIndexHandler<unsigned char>(0, 2),
                                                     FileChannelMap);

  // Retrieve the actual noise information.
  std::vector<double> PackedArray;
  for(size_t f = 0; f < 1024; f++) {
    NoiseMatrix& mat = noise.GetMatrixForIndex(f);
    const NoiseCorrelations::NoiseBlockIndexT& NoiseBlockIndex = noise.GetNoiseBlockIndex();
    assert(NoiseBlockIndex.MaxIndex() % 2 == 0 and NoiseBlockIndex.MaxIndex() > 0);

    // Create the name for this dataset.
    std::ostringstream strstream;
    strstream << "/noise_corr_" << std::setfill('0') << std::setw(4) << f;
    std::string dataset_name = strstream.str();

    // Get data from file.
    hid_t datasetID = H5Dopen2(fileID, dataset_name.c_str(), H5P_DEFAULT);
    hid_t dataspaceID = H5Dget_space(datasetID);
    size_t PackedSize = ExpectedPackedSize(f, FileBlockIndex.MaxIndex());
    assert(PackedSize == H5Sget_simple_extent_npoints(dataspaceID));
    PackedArray.resize(PackedSize);
    H5Dread(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            reinterpret_cast<void*>(&PackedArray[0]));

    // Fill noise matrix from packed array.
    for(size_t i = 0; i < NoiseBlockIndex.MaxIndex(); i++) {
      size_t rowInFile = FileBlockIndex.IndexForKey(NoiseBlockIndex.KeyForIndex(i));
      for(size_t j = i; j < NoiseBlockIndex.MaxIndex(); j++) {
        size_t colInFile = FileBlockIndex.IndexForKey(NoiseBlockIndex.KeyForIndex(j));

        std::pair<size_t, bool> packed_loc = PackedArrayIndexFor(f, rowInFile, colInFile, FileBlockIndex);
        if(packed_loc.first == size_t(-1)) mat.GetCorrByIndex(i, j) = 0;

        else {
          assert(packed_loc.first < PackedArray.size());
          if(packed_loc.second) mat.GetCorrByIndex(i, j) = -PackedArray[packed_loc.first];
          else mat.GetCorrByIndex(i, j) = PackedArray[packed_loc.first];
        }
      }
    }

    // Release HDF5 file resources.
    H5Sclose(dataspaceID);
    H5Dclose(datasetID);
  }

  assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
  H5Fclose(fileID);
}

} // namespace NoiseCorrelationIO
#endif
