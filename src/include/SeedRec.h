#ifndef SEEDREC_H
#define SEEDREC_H

//#include "InfoLists.h"
#include "SacRec.h"
#include <iostream>
#include <sstream>
#include <cstring>
#include <memory>

class SeedRec {

public:

   /* constructor: read in parameters through the input param file */
   SeedRec( const std::string& seedname, std::ostream& reportin );
   SeedRec( const std::string seedname = "", const std::string rdsexe = "",
				std::ostream& reportin = std::cerr );
   /* copy constructor */
   SeedRec( const SeedRec& SRin );
   /* move constructor */
   SeedRec( SeedRec&& SRin );
   /* assignment operator */
   SeedRec& operator= ( const SeedRec& SRin );
   /* move assignment operator */
	SeedRec& operator= ( SeedRec&& SRin );
   /* destructor */
   ~SeedRec();

   /* extracting sac from the seed file. */
	/*
	bool ExtractSac( const DailyInfo di, float& gapfrac, SacRec& sacout ) {
		return ExtractSac( di.staname, di.chname, di.sps, gapfrac, di.rec_outname, di.resp_outname, sacout );
	}
	*/
	bool ExtractSac( const std::string staname, const std::string chname, const int sps,
						  const std::string rec_outname, const std::string resp_outname, 
						  float& gapfrac, SacRec& sacout );

private:
   struct SRimpl;
   std::unique_ptr<SRimpl> pimpl;
	std::ostream* report = &(std::cerr);
};


#endif
