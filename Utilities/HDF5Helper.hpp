#ifndef HDF5Helper_hpp
#define HDF5Helper_hpp

/*
A common place to put hdf5 helper stuff.
For now, I install a callback routine which prints an error message and exits
whenever any hdf5 function experiences a failure;
callback is easier than checking the return value of all hdf5 functions.

There is a C++ interface to hdf5, but it is expected to undergo significant revision,
so I don't think it's worth adopting it yet.  If it were usable and offered exceptions,
that certainly would be better.
*/

#include "Utilities/IndexHandler.hpp"
#include "hdf5.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>
#include <string>

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

namespace HDF5Helper {

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

// We read/write APD channel maps a lot.
void WriteMapAsAttribute(const MapIndexHandler<unsigned char>& index, hid_t locID, std::string name) {
  assert(index.MaxIndex() > 0);
  const hsize_t size = index.MaxIndex();
  hid_t vectorID = H5Screate_simple(1, &size, NULL);
  hid_t attID = H5Acreate2(locID, name.c_str(), H5T_STD_U8LE, vectorID, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(&index.Buffer()));
  H5Aclose(attID);
  H5Sclose(vectorID);
}
MapIndexHandler<unsigned char> ReadMapFromAttribute(hid_t locID, std::string name) {
  MapIndexHandler<unsigned char> ret;
  hid_t attID = H5Aopen(locID, name.c_str(), H5P_DEFAULT);
  hid_t spaceID = H5Aget_space(attID);
  std::vector<unsigned char> tmpBuffer;
  tmpBuffer.resize(H5Sget_simple_extent_npoints(spaceID));
  H5Aread(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(&tmpBuffer[0]));
  for(size_t i = 0; i < tmpBuffer.size(); i++) ret.InsertKey(tmpBuffer[i]);
  assert(ret.MaxIndex() == tmpBuffer.size());
  H5Sclose(spaceID);
  H5Aclose(attID);
  return ret;
}

} // namespace HDF5Helper
#endif
