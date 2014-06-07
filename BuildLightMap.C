/*
Build the lightmap.  We take as input a pre-built sqlite database generated
by the BuildThEventDatabase.py script.

Note: if the sqlite queries are slow, we would probably benefit from compiling sqlite3 with the
SQLITE_ENABLE_STAT4 macro enabled before running ANALYZE.  However, we do the analysis from python's
sqlite3 module right now, so we don't have access to compilation options -- but there
are workarounds if necessary.
*/

#include "LightMap/LightMap.hpp"
#include "LightMap/LightMapIO.hpp"
#include "Utilities/IndexHandler.hpp"
#include "External/sqlite3.c"
#include <sstream>
#include <iostream>
#include <cstdlib>


int main()
{
  // Build an APD index.
  MapIndexHandler<unsigned char> APDIndex;
  for(unsigned char i = 152; i < 226; i++) {
    // Skip APD channels which have never been active.
    if(i == 178 or i == 191 or i == 205) continue;
    APDIndex.InsertKey(i);
  }

  // Build a lightmap object in memory.
  // This is where we specify the lightmap binning and choice of run ranges for gain.
  LightMap::PositionFunc posFunc;
  posFunc.SetAPDIndex(APDIndex);
  posFunc.SetBinning(-210, 15, 28, -210, 15, 28, -200, 10, 40); // 1.5cm in X/Y, 1cm in Z.
  LightMap::GainFunc gainFunc;
  gainFunc.InsertGain(APDIndex, 2401, 2423); // Old 09-28-11 APD biases
  gainFunc.InsertGain(APDIndex, 2424, 2690); // FEC voltage adjustment
  gainFunc.InsertGain(APDIndex, 2691, 2852); // Ebox 1 fan installed
  gainFunc.InsertGain(APDIndex, 2853, 2891); // Ebox 2 fan installed
  gainFunc.InsertGain(APDIndex, 2892, 3117); // Power outage here.
  gainFunc.InsertGain(APDIndex, 3118, 3326); // APD board swap
  gainFunc.InsertGain(APDIndex, 3327, 3700); // Utility power swap
  gainFunc.InsertGain(APDIndex, 3701, 3949);
  gainFunc.InsertGain(APDIndex, 3950, 4140); // Ralph's diode box installed
  gainFunc.InsertGain(APDIndex, 4141, 4579);
  gainFunc.InsertGain(APDIndex, 4580, 4779);
  gainFunc.InsertGain(APDIndex, 4780, 5197); // LC filters removed from FECs
  gainFunc.InsertGain(APDIndex, 5198, 5590); // Toward end of 5590 CR temps are elevated; APD gain impacted
  gainFunc.InsertGain(APDIndex, 5591, 5892); // Run2c ends.

  // Retrieve the indices of the x, y, and z coordinates.
  const IntervalIndexHandler& xIndex = posFunc.PosIndex().MajorIndex().MajorIndex();
  const IntervalIndexHandler& yIndex = posFunc.PosIndex().MajorIndex().MinorIndex();
  const IntervalIndexHandler& zIndex = posFunc.PosIndex().MinorIndex();

  // Initialize the gain function to be identically one.
  for(size_t i = 0; i < gainFunc.NumSnapshots(); i++) {
    LightMap::GainSnapshot& gain = gainFunc.GainAtIndex(i);
    for(size_t j = 0; j < gain.APDIndex().MaxIndex(); j++) gain.GetValAt(j) = 1;
  }

  // Create a connection to the database of thorium events.
  int ret;
  sqlite3* connection;
  ret = sqlite3_open_v2("Tmp/ThoriumLightmapEvents.db",
                        &connection,
                        SQLITE_OPEN_READONLY,
                        NULL);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to open database of events." << std::endl;
    exit(1);
  }

  // Build a comma-separated list of APD columns to retrieve in each case.
  std::ostringstream QueryColumns;
  QueryColumns << "apd_" << APDIndex.KeyForIndex(0) << "_magnitude";
  for(size_t i = 1; i < APDIndex.MaxIndex(); i++) {
    QueryColumns << ", apd_" << APDIndex.KeyForIndex(i) << "_magnitude";
  }
  QueryColumns << " FROM events WHERE ";

  // Prepare statements to retrieve by position or by run range.
  sqlite3_stmt* prep_stmt_pos, prep_stmt_run;
  ret = sqlite3_prepare_v2(connection,
                           std::string("SELECT runNo, " + QueryColumns.str() +
                                       "? <= xpos AND xpos < ? AND " +
                                       "? <= ypos AND ypos < ? AND " +
                                       "? <= zpos AND zpos < ?").c_str(),
                           -1,
                           &prep_stmt_pos,
                           NULL);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to prepare the position select statement." << std::endl;
    exit(1);
  }
  ret = sqlite3_prepare_v2(connection,
                           std::string("SELECT xpos, ypos, zpos, " + QueryColumns.str() +
                                       "? <= runNo AND runNo <= ?").c_str(),
                           -1,
                           &prep_stmt_run,
                           NULL);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to prepare the run select statement." << std::endl;
    exit(1);
  }

  // Iteratively converge on an appropriate lightmap -- three iterations seem to be enough.
  for(size_t iteration = 0; iteration < 3; iteration++) {

    std::cout << "Computing the position function, holding gain fixed." << std::endl;

    for(size_t posBin = 0; posBin < posFunc.PosIndex().MaxIndex(); posBin++) {

      // Get information about this position bin.
      size_t xyBin = posFunc.PosIndex().MajorIndexForIndex(posBin);
      size_t xBin = posFunc.PosIndex().MajorIndex().MajorIndexForIndex(xyBin);
      size_t yBin = posFunc.PosIndex().MajorIndex().MinorIndexForIndex(xyBin);
      size_t zBin = posFunc.PosIndex().MinorIndexForIndex(posBin);
      double xmin = xIndex.KeyForIndex(xBin);
      double xmax = xIndex.KeyForIndex(xBin+1);
      double ymin = yIndex.KeyForIndex(yBin);
      double ymax = yIndex.KeyForIndex(yBin+1);
      double zmin = zIndex.KeyForIndex(zBin);
      double zmax = zIndex.KeyForIndex(zBin+1);

      // Bind these parameters to the query.
      bool success = (
        sqlite3_bind_double(prep_stmt_pos, 1, xmin) == SQLITE_OK and
        sqlite3_bind_double(prep_stmt_pos, 2, xmax) == SQLITE_OK and
        sqlite3_bind_double(prep_stmt_pos, 3, ymin) == SQLITE_OK and
        sqlite3_bind_double(prep_stmt_pos, 4, ymax) == SQLITE_OK and
        sqlite3_bind_double(prep_stmt_pos, 5, zmin) == SQLITE_OK and
        sqlite3_bind_double(prep_stmt_pos, 6, zmax) == SQLITE_OK);
      if(not success) {
        std::cerr << "Failed to bind values to a select statement." << std::endl;
        exit(1);
      }

      // Retrieve all events which fall in this bin.
      // For each APD, we wish to compute the weighted sum of magnitudes
      // and the sum of square weights.
      LightMap::FuncVsAPD WeightedSumMagnitudes;
      LightMap::FuncVsAPD SumSquareWeights;
      for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
        WeightedSumMagnitudes[i] = 0;
        SumSquareWeights[i] = 0;
      }
      while((ret = sqlite3_step(prep_stmt_pos)) == SQLITE_ROW) {
        // Handle a row.
        int runNo = sqlite3_column_int(prep_stmt_pos, 0);
        const LightMap::GainSnapshot& gain = gainFunc.GainForRun(runNo);
        for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
          // Handle one APD for that row.
          double weight = gain.GetValAt(i);
          double apd_magnitude = sqlite3_column_double(prep_stmt_pos, 1 + i);
          WeightedSumMagnitudes[i] += weight*apd_magnitude;
          SumSquareWeights[i] += weight*weight;
        }
      }
      if(ret != SQLITE_DONE) {
        std::cerr << "Query returned failure." << std::endl;
        exit(1);
      }
      ret = sqlite3_reset(prep_stmt_pos);
      if(ret != SQLITE_OK) {
        std::cerr << "Failed to reset a statement." << std::endl;
        exit(1);
      }

      // Update the lightmap values for this bin, overwriting the old ones.
      for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
        if(WeightedSumMagnitudes[i] <= 0) posFunc.GetValAt(posBin, i) = 0;
        else posFunc.GetValAt(posBin, i) = WeightedSumMagnitudes[i]/SumSquareWeights[i];
      }
    }

    std::cout << "Done computing the position function." << std::endl;
    std::cout << "Computing the gain function, holding the position function fixed." << std::endl;

    for(size_t gainBin = 0; gainBin < gainFunc.NumSnapshots(); gainBin++) {
      LightMap::GainSnapshot& gain = gainFunc.GainAtIndex(gainBin);

      // Bind the run range parameters to the query.
      bool success = (
        sqlite3_bind_int(prep_stmt_run, 1, gain.FirstRun()) == SQLITE_OK and
        sqlite3_bind_int(prep_stmt_pos, 2, gain.LastRun()) == SQLITE_OK);
      if(not success) {
        std::cerr << "Failed to bind values to a select statement." << std::endl;
        exit(1);
      }

      // Retrieve all events which fall in this bin.
      // For each APD, we wish to compute the weighted sum of magnitudes
      // and the sum of square weights.
      LightMap::FuncVsAPD WeightedSumMagnitudes;
      LightMap::FuncVsAPD SumSquareWeights;
      for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
        WeightedSumMagnitudes[i] = 0;
        SumSquareWeights[i] = 0;
      }
      while((ret = sqlite3_step(prep_stmt_run)) == SQLITE_ROW) {
        // Handle a row.
        double xpos = sqlite3_column_double(prep_stmt_run, 0);
        double ypos = sqlite3_column_double(prep_stmt_run, 1);
        double zpos = sqlite3_column_double(prep_stmt_run, 2);

        size_t posBin = posFunc.PosIndex().IndexForKey(std::make_pair(std::make_pair(xpos, ypos), zpos));
        if(posBin >= posFunc.PosIndex().MaxIndex()) continue; // Kills events out of boundary for lightmap.
        const LightMap::FuncVsAPD& posFuncAtBin = posFunc.GetAllValsAt(posBin);

        for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
          // Handle one APD for that row.
          double weight = posFuncAtBin[i];
          double apd_magnitude = sqlite3_column_double(prep_stmt_run, 3 + i);
          WeightedSumMagnitudes[i] += weight*apd_magnitude;
          SumSquareWeights[i] += weight*weight;
        }
      }
      if(ret != SQLITE_DONE) {
        std::cerr << "Query returned failure." << std::endl;
        exit(1);
      }
      ret = sqlite3_reset(prep_stmt_run);
      if(ret != SQLITE_OK) {
        std::cerr << "Failed to reset a statement." << std::endl;
        exit(1);
      }

      // Update the gain snapshot, overwriting the old ones.
      for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
        if(WeightedSumMagnitudes[i] <= 0) gain.GetValAt(i) = 0;
        else gain.GetValAt(i) = WeightedSumMagnitudes[i]/SumSquareWeights[i];
      }
    }
    std::cout << "Done computing the gain function." << std::endl;
  } // End loop over iterations.

  // Clean up sqlite resources.
  ret = sqlite3_finalize(prep_stmt_pos);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to finalize an sqlite statement." << std::endl;
    std::cerr << "We'll finish up because we (probably) can, but this should be fixed." << std::endl;
  }
  ret = sqlite3_finalize(prep_stmt_run);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to finalize an sqlite statement." << std::endl;
    std::cerr << "We'll finish up because we (probably) can, but this should be fixed." << std::endl;
  }
  ret = sqlite3_close_v2(connection);
  if(ret != SQLITE_OK) {
    std::cerr << "Failed to close the sqlite connection." << std::endl;
    std::cerr << "We'll finish up because we (probably) can, but this should be fixed." << std::endl;
  }

  // We've generated a good (hopefully) lightmap -- now write it to file.
  std::cout << "Finished generating lightmap -- write to Tmp/LightMap.hdf5." << std::endl;
  LightMapIO::WriteLightMap("Tmp/LightMap.hdf5", posFunc, gainFunc);
  std::cout << "Done writing the lightmap to file.  Have a nice day." << std::endl;
}
