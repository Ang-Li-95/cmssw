#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"

#include "RecoLocalFastTime/FTLClusterizer/interface/MTDCPEFromSiPMTime.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"
// MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"


#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------
//  A constructor run for generic and templates
//
//-----------------------------------------------------------------------------
MTDCPEFromSiPMTime::MTDCPEFromSiPMTime(edm::ParameterSet const & conf,
		       		       const MTDGeometry& geom)
  : MTDCPEBase(conf, geom)
{
   //fillDetParams();
}


//-----------------------------------------------------------------------------
//  Fill all variables which are constant for an event (geometry)
//-----------------------------------------------------------------------------
/*
void MTDCPEFromSiPMTime::fillDetParams()
{
  auto const & dus = geom_.detUnits();
  unsigned detectors = dus.size();
  m_DetParams.resize(detectors);
  LogDebug("MTDCPEBase::fillDetParams():") <<"caching "<<detectors<<"MTD detectors"<<endl;
  for (unsigned i=0; i!=detectors;++i) 
    {
      auto & p=m_DetParams[i];
      p.theDet = dynamic_cast<const MTDGeomDetUnit*>(dus[i]);
      assert(p.theDet);
      
      p.theOrigin = p.theDet->surface().toLocal(GlobalPoint(0,0,0));
      
      //--- p.theDet->type() returns a GeomDetType, which implements subDetector()
      p.thePart = p.theDet->type().subDetector();
            
      //--- bounds() is implemented in BoundSurface itself.
      p.theThickness = p.theDet->surface().bounds().thickness();
      
      // Cache the det id for templates and generic erros
      p.theTopol = &(static_cast<const ProxyMTDTopology&>(p.theDet->topology()));
      assert(p.theTopol);
      p.theRecTopol = &(static_cast<const RectangularMTDTopology&>(p.theTopol->specificTopology()));    
      assert(p.theRecTopol);
      
      //--- The geometrical description of one module/plaquette
      std::pair<float,float> pitchxy = p.theRecTopol->pitch();
      p.thePitchX = pitchxy.first;	     // pitch along x
      p.thePitchY = pitchxy.second;	     // pitch along y
            
      LogDebug("MTDCPEBase::fillDetParams()") << "***** MTD LAYOUT *****"
					      << " thePart = " << p.thePart
					      << " theThickness = " << p.theThickness
					      << " thePitchX  = " << p.thePitchX
					      << " thePitchY  = " << p.thePitchY;
    }
}

//-----------------------------------------------------------------------------
//  One function to cache the variables common for one DetUnit.
//-----------------------------------------------------------------------------
void
MTDCPEBase::setTheClu( DetParam const & dp, ClusterParam & cp ) const
{   
}

//------------------------------------------------------------------------
MTDCPEBase::DetParam const & MTDCPEBase::detParam(const GeomDetUnit & det) const 
{
  return m_DetParams.at(det.index());
}*/
/*
LocalPoint
MTDCPEFromSiPMTime::localPosition(DetParam const & dp, ClusterParam & cp) const
{
  //remember measurement point is row(col)+0.5f
  MeasurementPoint pos(cp.theCluster->x(),cp.theCluster->y());
  return dp.theTopol->localPosition(pos);
}
*/
LocalError
MTDCPEFromSiPMTime::localError(DetParam const & dp,  ClusterParam & cp) const
{
  constexpr double one_over_twelve = 1./12.;
  //constexpr double one_over_twelve = 10;
  MeasurementPoint pos(cp.theCluster->x(),cp.theCluster->y());
  MeasurementError simpleRect(one_over_twelve,0,one_over_twelve);
  if (GeomDetEnumerators::isBarrel(dp.thePart) ){
    MeasurementError simpleRect_new(one_over_twelve,0,0.0164);
    simpleRect = simpleRect_new;
  }
  /*
  else {
  if (GeomDetEnumerators::isEndcap(dp.thePart) ){
    MeasurementError simpleRect(one_over_twelve,0,one_over_twelve);
  }*/
  return dp.theTopol->localError(pos,simpleRect);


}
/*
MTDCPEFromSiPMTime::TimeValue
MTDCPEFromSiPMTime::clusterTime(DetParam const & dp, ClusterParam & cp) const
{
  return cp.theCluster->time();
}


MTDCPEFromSiPMTime::TimeValueError
MTDCPEFromSiPMTime::clusterTimeError(DetParam const & dp, ClusterParam & cp) const
{
  return cp.theCluster->timeError();
}*/
