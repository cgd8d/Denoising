#ifndef LightMapIO_hpp
#define LightMapIO_hpp

#include "LightMap/LightMap.hpp"
#include "Utilities/HDF5Helper.hpp"
#include "hdf5.h"
#include <string>


namespace LightMapIO {

/*
Write the position and gain functions to a single hdf5 file.
If the file already exists, fail and exit.
*/
void WriteLightMap(std::string filename,
                   const LightMap::PositionFunc& posFunc,
                   const LightMap::GainFunc& gainFunc)
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

  // Write MAX_APDS as well -- this should never change, but we will verify it.
  {
    hid_t scalarID = H5Screate(H5S_SCALAR);
    hid_t attID = H5Acreate2(fileID, "MAX_APDS", H5T_STD_U8LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<void*>(&LightMap::MAX_APDS));
    H5Aclose(attID);
    H5Sclose(scalarID);
  }

  // Write the position function in the root group.
  // Also attach the APD index.
  {
    hsize_t SizeOfArray = hsize_t(LightMap::MaxAPDS)*posFunc.PosIndex().MaxIndex();
    hid_t vectorID = H5Screate_simple(1, &SizeOfArray, NULL);
    hid_t datasetID = H5Dcreate2(fileID, "posFunc", H5T_IEEE_F64LE, vectorID,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             reinterpret_cast<void*>(&posFunc.GetValAt(0, 0));
    HDF5Helper::WriteMapAsAttribute(posFunc.APDIndex(), datasetID, "APDindex");
    H5Dclose(datasetID);
    H5Sclose(vectorID);
  }

  // Write gain functions.
  {
    hid_t grpID = H5Gcreate2(fileID, "/gainFunc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(size_t i = 0; i < gainFunc.NumSnapshots(); i++) {
      const GainSnapshot& gain = gainFunc.GainAtIndex(i);

      // Produce the appropriate name.
      std::ostringstream namestr;
      namestr << "gain_" << gain.FirstRun() << "_" << gain.LastRun();
      std::string name = namestr.str();

      // Write stuff to file.
      hsize_t SizeOfArray = gain.APDIndex().MaxIndex();
      hid_t vectorID = H5Screate_simple(1, &SizeOfArray, NULL);
      hid_t datasetID = H5Dcreate2(grpID, name.str(), H5T_IEEE_F64LE, vectorID,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               reinterpret_cast<void*>(&gain.GetValAt(0));
      HDF5Helper::WriteMapAsAttribute(gain.APDIndex(), datasetID, "APDindex");
      int firstRun = gain.FirstRun();
      int lastRun = gain.LastRun();
      hid_t scalarID = H5Screate(H5S_SCALAR);
      hid_t firstID = H5Acreate2(datasetID, "first_run", H5T_STD_I32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(firstID, H5T_NATIVE_INT, reinterpret_cast<void*>(&firstRun));
      hid_t lastID = H5Acreate2(datasetID, "last_run", H5T_STD_I32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(lastID, H5T_NATIVE_INT, reinterpret_cast<void*>(&lastRun));
      H5Aclose(lastID);
      H5Aclose(firstID);
      H5Sclose(scalarID);
      H5Dclose(datasetID);
      H5Sclose(vectorID);
    }
    H5Gclose(grpID);
  }

  assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 0); // Make sure we didn't leak any objects.
  H5Fclose(fileID);
}

} // namespace LightMapIO
#endif
