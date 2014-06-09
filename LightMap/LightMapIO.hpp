#ifndef LightMapIO_hpp
#define LightMapIO_hpp

/*
Perform IO on the various lightmap objects, reading from and writing to hdf5 files.

There are many ancillary functions (and even one object) here; the user will generally need:

namespace LightMapIO {
void WriteLightMap(std::string filename, const LightMap::PositionFunc& posFunc, const LightMap::GainFunc& gainFunc);
LightMap::GainSnapshot ReadLightmapAtRun(std::string filename, int run, LightMap::PositionFunc& posFunc);
LightMap::GainFunc ReadLightmap(std::string filename, LightMap::PositionFunc& posFunc);
class GainFuncNotKnown;
}

These functions write a lightmap to an hdf5 file; read it out again (either reading all gain values,
or only the gain corresponding to a particular run); and, in the absence of a gain snapshot
matching a provided run number, raise a GainFuncNotKnown exception.
*/

#include "LightMap/LightMap.hpp"
#include "Utilities/HDF5Helper.hpp"
#include "Utilities/IndexHandler.hpp"
#include "hdf5.h"
#include <string>
#include <sstream>
#include <stdexcept>

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
    const unsigned char version = 1;
    hid_t scalarID = H5Screate(H5S_SCALAR);
    hid_t attID = H5Acreate2(fileID, "version", H5T_STD_U8LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<const void*>(&version));
    H5Aclose(attID);
    H5Sclose(scalarID);
  }

  // Write MAX_APDS as well -- this should never change, but we will verify it.
  {
    hid_t scalarID = H5Screate(H5S_SCALAR);
    hid_t attID = H5Acreate2(fileID, "MAX_APDS", H5T_STD_U8LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_UCHAR, reinterpret_cast<const void*>(&LightMap::MAX_APDS));
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
             reinterpret_cast<const void*>(&posFunc.GetValAt(0, 0)));
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
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xmin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "xstep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xstep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "xnum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&xnum));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ymin", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ymin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ystep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&ystep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "ynum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&ynum));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "zmin", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&xmin));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "zstep", H5T_IEEE_F64LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_DOUBLE, reinterpret_cast<void*>(&zstep));
    H5Aclose(attID);
    attID = H5Acreate2(datasetID, "znum", H5T_STD_U32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attID, H5T_NATIVE_HSIZE, reinterpret_cast<void*>(&znum));
    H5Aclose(attID);

    H5Sclose(scalarID);
    H5Dclose(datasetID);
    H5Sclose(vectorID);
  }

  // Write gain functions.
  {
    hid_t grpID = H5Gcreate2(fileID, "/gainFunc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(size_t i = 0; i < gainFunc.NumSnapshots(); i++) {
      const LightMap::GainSnapshot& gain = gainFunc.GainAtIndex(i);

      // Produce the appropriate name.
      std::ostringstream namestr;
      namestr << "gain_" << gain.FirstRun() << "_" << gain.LastRun();
      std::string name = namestr.str();

      // Write stuff to file.
      hsize_t SizeOfArray = gain.APDIndex().MaxIndex();
      hid_t vectorID = H5Screate_simple(1, &SizeOfArray, NULL);
      hid_t datasetID = H5Dcreate2(grpID, name.c_str(), H5T_IEEE_F64LE, vectorID,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               reinterpret_cast<const void*>(&gain.GetValAt(0)));
      HDF5Helper::WriteMapAsAttribute(gain.APDIndex(), datasetID, "APDindex");
      const int firstRun = gain.FirstRun();
      const int lastRun = gain.LastRun();
      hid_t scalarID = H5Screate(H5S_SCALAR);
      hid_t firstID = H5Acreate2(datasetID, "first_run", H5T_STD_I32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(firstID, H5T_NATIVE_INT, reinterpret_cast<const void*>(&firstRun));
      hid_t lastID = H5Acreate2(datasetID, "last_run", H5T_STD_I32LE, scalarID, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(lastID, H5T_NATIVE_INT, reinterpret_cast<const void*>(&lastRun));
      H5Aclose(lastID);
      H5Aclose(firstID);
      H5Sclose(scalarID);
      H5Dclose(datasetID);
      H5Sclose(vectorID);
    }
    H5Gclose(grpID);
  }

  assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
  H5Fclose(fileID);
}

/*
Read in the position function from an already-opened hdf5 file.
Assume global attributes have already been checked.
*/
void ReadPosFunc(hid_t fileID, LightMap::PositionFunc& posFunc) {
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

/*
Read out a gain snapshot.  This assumes you've already found the HDF5 dataset corresponding
to the snapshot you want.  It is not meant to be called externally.
Note that the cost of copying a gain snapshot is small, and there aren't many,
so even if the compiler doesn't elide the copy it's OK.
*/
LightMap::GainSnapshot ReadGainSnapshotAt(hid_t datasetID) {

  // First get the attributes.
  LightMap::APDIndexT apd_index = HDF5Helper::ReadMapFromAttribute(datasetID, "APDindex");
  hid_t attID;
  int first_run, last_run;
  attID = H5Aopen(datasetID, "first_run", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_INT, reinterpret_cast<void*>(&first_run));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "last_run", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_INT, reinterpret_cast<void*>(&last_run));
  H5Aclose(attID);
  assert(first_run <= last_run);
  assert(apd_index.MaxIndex() <= LightMap::MAX_APDS);

  // Now make the gain snapshot and read in the data.
  LightMap::GainSnapshot gain(apd_index, first_run, last_run);
  hid_t dataspaceID = H5Dget_space(datasetID);
  assert(apd_index.MaxIndex() == H5Sget_simple_extent_npoints(dataspaceID));
  H5Dread(datasetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          reinterpret_cast<void*>(&gain.GetValAt(0)));
  H5Sclose(dataspaceID);

  return gain;
}

/*
Find the snapshot corresponding to a particular run, and return it.
This requires iterating through the group.
We require the user to hand us a file id from an already-open HDF5 file.
If the run has no matching snapshot, raise an exception -- we'll probably
choose to use a default gain map, with reservations.
*/
class GainFuncNotKnown : public std::runtime_error
{
 public:
  // Identify errors due to a file not containing an appropriate gain function.
  // We'll generally recover from this circumstance by generating a default
  // gain function, but we'll need to know it happened well upstream of here.
  GainFuncNotKnown(const std::string& what) : std::runtime_error(what) { }
};

herr_t TestSnapshotContainsRun(hid_t g_id, const char *name,
                               const H5L_info_t * /* unused */, void *op_data) {
  int* runNo = reinterpret_cast<int*>(op_data); // We're handed the run number we seek.
  hid_t datasetID = H5Dopen2(g_id, name, H5P_DEFAULT);
  hid_t attID;
  int first_run, last_run;
  attID = H5Aopen(datasetID, "first_run", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_INT, reinterpret_cast<void*>(&first_run));
  H5Aclose(attID);
  attID = H5Aopen(datasetID, "last_run", H5P_DEFAULT);
  H5Aread(attID, H5T_NATIVE_INT, reinterpret_cast<void*>(&last_run));
  H5Aclose(attID);
  assert(first_run <= last_run);

  // If we find the right object, don't close it!
  // It will have a positive ID; return it to finish iteration.
  // Note the conversion from hid_t to herr_t and back -- perhaps not totally kosher.
  if(first_run <= *runNo and *runNo <= last_run) return datasetID;

  // Otherwise, close the dataset and continue.
  H5Dclose(datasetID);
  return 0;
}

LightMap::GainSnapshot ReadGainSnapshotForRun(int run, hid_t fileID) {
  hid_t grpID = H5Gopen2(fileID, "/gainFunc", H5P_DEFAULT);

  // Iterate through the group, looking for the right gain snapshot.
  hid_t ret = (hid_t)H5Literate(grpID, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
                                 &TestSnapshotContainsRun, reinterpret_cast<void*>(&run));
  if(ret == 0) {
    H5Gclose(grpID); // Don't leak resources.
    throw GainFuncNotKnown("This run does not have a known gain function.");
  }

  LightMap::GainSnapshot gain = ReadGainSnapshotAt(ret);
  H5Dclose(ret);
  H5Gclose(grpID);
  return gain;
}

/*
Read all of the lightmap information for a particular run from a file.
The PosFunc will be filled; we also return a gain snapshot.
If the gain snapshot does not exist in the file, we fill the posFunc and
throw an exception of type GainFuncNotKnown.
*/
LightMap::GainSnapshot ReadLightmapAtRun(std::string filename, int run, LightMap::PositionFunc& posFunc) {
  hid_t fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the position function.
  ReadPosFunc(fileID, posFunc);

  // Read the gain snapshot.
  try {
    LightMap::GainSnapshot gain = ReadGainSnapshotForRun(run, fileID);
    assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
    H5Fclose(fileID);
    return gain;
  }
  catch (GainFuncNotKnown& exc) {
    assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
    H5Fclose(fileID);
    throw; // Re-throw the exception -- just needed to clean up the file.
  }
}

/*
Read the whole lightmap, including *all* gain snapshots.
Fill the posFunc which is passed in; return a GainFunc as well.
*/
herr_t AddSnapshotToGainFunc(hid_t g_id, const char *name,
                             const H5L_info_t * /* unused */, void *op_data) {
  // Adaptor from HDF5 callback signature.
  LightMap::GainFunc* gainFunc = reinterpret_cast<LightMap::GainFunc*>(op_data);
  hid_t datasetID = H5Dopen2(g_id, name, H5P_DEFAULT);
  gainFunc->InsertGain(ReadGainSnapshotAt(datasetID));
  H5Dclose(datasetID);
  return 0;
}

LightMap::GainFunc ReadLightmap(std::string filename, LightMap::PositionFunc& posFunc) {
  hid_t fileID = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  // Read the position function.
  ReadPosFunc(fileID, posFunc);

  // Fill a GainFunc as well.
  LightMap::GainFunc gainFunc;
  hid_t grpID = H5Gopen2(fileID, "/gainFunc", H5P_DEFAULT);
  H5Literate(grpID, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
             &AddSnapshotToGainFunc, reinterpret_cast<void*>(&gainFunc));

  // Clean up and return.
  H5Gclose(grpID);
  assert(H5Fget_obj_count(fileID, H5F_OBJ_ALL) == 1); // The file should be the only object left.
  H5Fclose(fileID);
  return gainFunc;
}

} // namespace LightMapIO
#endif
