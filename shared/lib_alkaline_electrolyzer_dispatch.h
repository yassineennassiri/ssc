#ifndef _LIB_ALKALINE_ELECTROLYZER_DISPATCH_
#DEFINE _LIB_ALKALINE_ELECTROLYZER_DISPATCH_

#include "lib_alkaline_electrolyzer.h"

class AlkalineElectrolyzerDispatch
{
public:
        /// Default constructor
        AlkalineElectrolyzerDispatch() { /* Nothing to do */ };
        
        /// Construct with argument
        AlkalineElectrolyzerDispatch();
        
        /// Destructor
        ~AlkalineElectrolyzer();
        
        /// Run dispatch for single step
        void runSingleTimeStep();
        
        /// Update dispatch units (for testing)
        void setDispatchOption();
        
        /// Update the fixed percentage to dispatch
        void setFixedDischargePercentage();
        
        /// Update dispatch units (for testing)
        void seManualDispacthUnits();
        
        /// Get the total Hydrogen Produced kg
