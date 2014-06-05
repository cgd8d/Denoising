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
  // Also attach the APD index and position indices.
  {
    hsize_t SizeOfArray = hsize_t(LightMap::MAX_APDS)*posFunc.PosIndex().MaxIndex();
    hid_t vectorID = H5Screate_simple(1, &SizeOfArray, NULL);
    hid_t datasetID = H5Dcreate2(fileID, "posFunc", H5T_IEEE_F64LE, vectorID,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             reinterpret_cast<void*>(&posFunc.GetValAt(0, 0));
    HDF5Helper::WriteMapAsAttribute(posFunc.APDIndex(), datasetID, "APDindex");

    const LightMap::PositionFunc::PosIndexT& pos_index = posFunc.PosIndex();
    double xmin = pos_index.MajorIndex().MajorIndex().Start();
    double xstep = pos_index.MajorIndex().MajorIndex().StepSize();
    size_t xnum = pos_index.MajorIndex().MajorIndex().MaxIndex();
    double ymin = pos_index.MajorIndex().MinorIndex().Start();
    double ystep = pos_index.MajorIndex().MinorIndex().StepSize();
    size_t ynum = pos_index.MajorIndex().MinorIndex().MaxIndex();
    double zmin = pos_index.MinorIndex().Start();
    double zstep = pos_index.MinorIndex().StepSize();
    size_t znum = pos_index.MinorIndex().MaxIndex();

    hid_t scalarID = H5Screate(H5S_SCALAR);
    hid_t attID;
    attID = H5Acreate2(datasetID, "xmin", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xmin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "xstep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xstep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "xnum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&xnum));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ymin", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ymin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ystep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ystep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ynum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&ynum));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "zmin", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xmin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "zstep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&zstep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "znum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(firstID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&znum));
    H5Aclose(attID);

    H5Sclose(scalarID);
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

/*
Read in the position function from an already-opened hdf5 file.
Assume global attributes have already been checked.
*/
void ReadPosFunc(hid_t fileID, LightMap::PosFunc& posFunc) {
  hid_t datasetID = H5Dopen2(fileID, "posFunc", H5P_DEFAULT);

  // Start by reading attributes.
  double xmin, xstep, ymin, ystep, zmin, zstep;
  size_t xnum, ynum, znum;
  hid_t attID;
  attID = H5Aopen(datasetID, "xmin", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xmin));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "xstep", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xstep));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "xnum", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&xnum));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "ymin", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ymin));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "ystep", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ystep));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "ynum", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&ynum));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "zmin", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&zmin));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "zstep", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&zstep));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "znum", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&znum));
  H5Aclose(attID);

  LightMap::APDIndexT apdIndex = HDF5Helper::ReadMapFromAttribute(datasetID, "APDindex");
  assert(apdIndex.MaxIndex() <= LightMap::MAX_APDS);

  posFunc.SetAPDIndex(apdIndex);
  posFunc.SetBinning(xmin, xstep, xnum, ymin, ystep, ynum, zmin, zstep, znum);

  // Now we know enough to read in the big chunk of data.
  hsize_t SizeOfArray = hsize_t(LightMap::MAX_APDS)*posFunc.PosIndex().MaxIndex();
  hid_t dataspaceID = H5Dget_space(datasetID);
  assert(SizeOfArray == H5Sget_simple_extent_npoints(dataspaceID));
  H5Dread(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          reinterpret_cast<void*>(&posFunc.GetValAt(0, 0)));

  // Close things to avoid a leak, and return.
  H5Sclose(dataspaceID);
  H5Dclose(datasetID);
}

} // namespace LightMapIO
#endif
