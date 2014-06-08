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
#include "External/sqlite3.h"
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
  sqlite3_stmt *prep_stmt_pos, *prep_stmt_run;
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

  // We know relative R(x) and S(t) -- fix them using [hard-coded] laser data.
  // The laser data used is from laser run 4540.
  // NOTE: by the time of run 4540, channel 163 was already disconnected, so the number here is made up.
  // NOTE: it would be great in the future to do a better correspondence between
  // absolute gain measurements from laser runs and relative gain measurements from thorium.
  for(size_t i = 0; i < APDIndex.MaxIndex(); i++) {
    double laser_gain;
    switch(APDIndex.KeyForIndex(i)) {
      case 152: laser_gain = 201.230438146; break;
      case 153: laser_gain = 178.750438779; break;
      case 154: laser_gain = 194.228589338; break;
      case 155: laser_gain = 183.33801615; break;
      case 156: laser_gain = 218.485999976; break;
      case 157: laser_gain = 222.139259152; break;
      case 158: laser_gain = 169.982559736; break;
      case 159: laser_gain = 140.385120552; break;
      case 160: laser_gain = 137.602725389; break;
      case 161: laser_gain = 197.78183714; break;
      case 162: laser_gain = 155.478773762; break;
      case 163: laser_gain = 200; break; // FICTITIOUS, but this channel was disconnected Feb 2012.  Better guess?
      case 164: laser_gain = 175.875067527; break;
      case 165: laser_gain = 160.014408865; break;
      case 166: laser_gain = 183.408055613; break;
      case 167: laser_gain = 189.600819126; break;
      case 168: laser_gain = 160.339214431; break;
      case 169: laser_gain = 168.547991045; break;
      case 170: laser_gain = 182.670039836; break;
      case 171: laser_gain = 205.567802982; break;
      case 172: laser_gain = 195.87450621; break;
      case 173: laser_gain = 224.956647122; break;
      case 174: laser_gain = 232.062359991; break;
      case 175: laser_gain = 241.822881767; break;
      case 176: laser_gain = 194.740435753; break;
      case 177: laser_gain = 189.867775084; break;
      // case 178: laser_gain = 0; // Bad channel, omitted.
      case 179: laser_gain = 206.755206938; break;
      case 180: laser_gain = 207.822617603; break;
      case 181: laser_gain = 207.501985741; break;
      case 182: laser_gain = 218.213137769; break;
      case 183: laser_gain = 234.369354843; break;
      case 184: laser_gain = 99.908111992; break;
      case 185: laser_gain = 238.381809313; break;
      case 186: laser_gain = 225.118270743; break;
      case 187: laser_gain = 199.078450518; break;
      case 188: laser_gain = 221.863823239; break;
      case 189: laser_gain = 177.032783679; break;
      case 190: laser_gain = 196.787332164; break;
      // case 191: laser_gain = 0; // Bad channel, omitted.
      case 192: laser_gain = 194.923448865; break;
      case 193: laser_gain = 197.027984846; break;
      case 194: laser_gain = 202.757086104; break;
      case 195: laser_gain = 194.432937658; break;
      case 196: laser_gain = 208.992809367; break;
      case 197: laser_gain = 224.762562055; break;
      case 198: laser_gain = 217.696006443; break;
      case 199: laser_gain = 222.380158829; break;
      case 200: laser_gain = 218.358804472; break;
      case 201: laser_gain = 209.573057132; break;
      case 202: laser_gain = 194.684536629; break;
      case 203: laser_gain = 182.543842783; break;
      case 204: laser_gain = 193.469930111; break;
      // case 205: laser_gain = 0; // Bad channel, omitted.
      case 206: laser_gain = 193.627191472; break;
      case 207: laser_gain = 196.073150574; break;
      case 208: laser_gain = 189.597962521; break;
      case 209: laser_gain = 198.824317108; break;
      case 210: laser_gain = 222.747770671; break;
      case 211: laser_gain = 216.928470825; break;
      case 212: laser_gain = 223.437239807; break;
      case 213: laser_gain = 224.316404923; break;
      case 214: laser_gain = 216.26783603; break;
      case 215: laser_gain = 209.612423384; break;
      case 216: laser_gain = 223.041660884; break;
      case 217: laser_gain = 202.642254512; break;
      case 218: laser_gain = 213.904993632; break;
      case 219: laser_gain = 221.988942321; break;
      case 220: laser_gain = 201.427174798; break;
      case 221: laser_gain = 196.689200146; break;
      case 222: laser_gain = 191.457656123; break;
      case 223: laser_gain = 186.183873541; break;
      case 224: laser_gain = 217.033080346; break;
      case 225: laser_gain = 205.858374653; break;
      default: laser_gain = 0; // Bad or non-existent channel.
    }

    double thorium_gain = gainFunc.GainForRun(4540).GetValAt(i);

    if(laser_gain == 0 or thorium_gain == 0) continue;

    double ratio_gains = laser_gain/thorium_gain;
    for(size_t j = 0; j < gainFunc.NumSnapshots(); j++) {
      gainFunc.GainAtIndex(j).GetValAt(i) *= ratio_gains;
    }
    for(size_t j = 0; j < posFunc.PosIndex().MaxIndex(); j++) {
      posFunc.GetValAt(j, i) /= ratio_gains; // To keep the product of R and S the same.
    }
  }

  // We've generated a good (hopefully) lightmap -- now write it to file.
  std::cout << "Finished generating lightmap -- write to Tmp/LightMap.hdf5." << std::endl;
  LightMapIO::WriteLightMap("Tmp/LightMap.hdf5", posFunc, gainFunc);
  std::cout << "Done writing the lightmap to file.  Have a nice day." << std::endl;
}
