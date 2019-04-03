#ifndef __LIB_ALKALINE_ELECTROLYZER__
#define __LIB_ALKALINE_ELECTROLYZER__

#include <map>
#include "lib_util.h"

/** 
* \ class AlkalineElectrolyzer
*
* \brief
*
*  The AlkalineElectrolyzer class provides the overall proporties and operation of the electrolyzer to model 
   the Alkaline electrolyzer technology in SAM.
*/

class AlkalineElectrolyzer
{
  public:
          /// Default AlkalineElectrolyzer constructor
          AlkalineElectrolyzer();
  
          /// Construct Alkaline Electrolyzer arguments
          AlkalineElectrolyzer('arguments') ;
  
  
          /// Default destructor
          ~AlkalineElectrolyzer();
  
  
          /// Copy Destructor 
          AlkalineElectrolyzer(const AlkalineElectrolyzer &alkalineElectrolyzer);
  
  
          /// Run for single time step
          void runSingleTimeStep (double power_kW);
   
          /// Return true if starting up but not fully running
          bool isStarting();
   
          /// Return true if operating 
          bool isRunning();
   
          /// Return true if shutting down
          bool isShuttingDown();
   
          /// Return true if totally shut down
          bool isShutDown();
   
          /// Return false if hour zero needs initial power from dispatch 
          bool isInitialized();
   
   
   
