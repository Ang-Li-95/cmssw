#ifndef RecoLocalFastTime_FTLClusterizer_MTDCPEFromSiPMTime_H
#define RecoLocalFastTime_FTLClusterizer_MTDCPEFromSiPMTime_H 1

//-----------------------------------------------------------------------------
// \class        MTDCPEFromSiPMTime
//-----------------------------------------------------------------------------

#include <utility>
#include <memory>
#include <vector>
#include "TMath.h"

#include "MTDClusterParameterEstimator.h"

#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

#include "RecoLocalFastTime/FTLClusterizer/interface/MTDCPEBase.h"

#include <unordered_map>

#include <iostream>


class MTDCPEFromSiPMTime : public MTDCPEBase
{
/*
public:
   struct DetParam
   {
      DetParam() {}
      const MTDGeomDetUnit * theDet;
      const ProxyMTDTopology * theTopol;
      const RectangularMTDTopology * theRecTopol;

      GeomDetType::SubDetector thePart;      
      Local3DPoint theOrigin;
      float theThickness;
      float thePitchX;
      float thePitchY;
   };
   
   struct ClusterParam
   {
      ClusterParam(const FTLCluster & cl) : theCluster(&cl) {}

      virtual ~ClusterParam() = default;

      const FTLCluster * theCluster;
   };
   */
   
public:
   MTDCPEFromSiPMTime(edm::ParameterSet const& conf, const MTDGeometry& geom);  
   
   /*
   inline ReturnType getParameters(const FTLCluster & cl,
                                   const GeomDetUnit    & det ) const override
   {
      
      DetParam const & dp = detParam(det);
      ClusterParam cp(cl);
      setTheClu( dp, cp );
      auto tuple = std::make_tuple(
				   localPosition(dp, cp),
				   localError(dp, cp),
				   clusterTime(dp, cp),
				   clusterTimeError(dp, cp)
				   );
      return tuple;
   }
   
   //--------------------------------------------------------------------------
   // In principle we could use the track too to get an angle if needed
   //--------------------------------------------------------------------------
   inline ReturnType getParameters(const FTLCluster & cl,
                                   const GeomDetUnit    & det,
                                   const LocalTrajectoryParameters & ltp ) const override
   {
     return getParameters(cl,det);
   }
   */
   
   
private:
   //--------------------------------------------------------------------------
   // This is where the action happens.
   //--------------------------------------------------------------------------
   virtual LocalPoint localPosition(DetParam const & dp, ClusterParam & cp) const override;
   virtual LocalError localError   (DetParam const & dp, ClusterParam & cp) const override;
   //virtual TimeValue  clusterTime(DetParam const & dp, ClusterParam & cp) const;
   //virtual TimeValueError  clusterTimeError(DetParam const & dp, ClusterParam & cp) const;
   
   //void fillDetParams();

/*   
protected:   
   //---------------------------------------------------------------------------
   //  Data members
   //---------------------------------------------------------------------------
   
   //--- Global quantities   
   const MTDGeometry & geom_;          // geometry
   
protected:
   
   void  setTheClu( DetParam const & dp, ClusterParam & cp ) const ;   
   
   //---------------------------------------------------------------------------
   //  Cluster-level services.
   //---------------------------------------------------------------------------
   
   DetParam const & detParam(const GeomDetUnit & det) const;
   
   using DetParams=std::vector<DetParam>;
   DetParams m_DetParams;
*/   
};

#endif


