// Our own includes
#include "RecoLocalFastTime/FTLClusterizer/interface/MTDArrayBuffer.h"
#include "RecoLocalFastTime/FTLClusterizer/interface/MTDThresholdClusterizer.h"

#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

// STL
#include <stack>
#include <vector>
#include <iostream>
#include <atomic>
#include <cmath>
using namespace std;

//#define DEBUG_ENABLED 
#ifdef DEBUG_ENABLED
#define DEBUG(x) do { std::cout << x << std::endl; } while (0)
#else
#define DEBUG(x)
#endif

//----------------------------------------------------------------------------
//! Constructor: 
//----------------------------------------------------------------------------
MTDThresholdClusterizer::MTDThresholdClusterizer
  (edm::ParameterSet const& conf) :
    // Get energy thresholds 
    theHitThreshold( conf.getParameter<double>("HitThreshold") ),
    theSeedThreshold( conf.getParameter<double>("SeedThreshold") ),
    theClusterThreshold( conf.getParameter<double>("ClusterThreshold") ),
    theTimeThreshold( conf.getParameter<double>("TimeThreshold") ),
    theNumOfRows(0), theNumOfCols(0), theCurrentId(0), 
    theBuffer(theNumOfRows, theNumOfCols ),
    bufferAlreadySet(false)
{
}
/////////////////////////////////////////////////////////////////////////////
MTDThresholdClusterizer::~MTDThresholdClusterizer() {}


// Configuration descriptions
void
MTDThresholdClusterizer::fillDescriptions(edm::ParameterSetDescription& desc) {
  desc.add<double>("HitThreshold", 0.);
  desc.add<double>("SeedThreshold", 0.);
  desc.add<double>("ClusterThreshold", 0.);
  desc.add<double>("TimeThreshold", 10.);
}

//----------------------------------------------------------------------------
//!  Prepare the Clusterizer to work on a particular DetUnit.  Re-init the
//!  size of the panel/plaquette (so update nrows and ncols), 
//----------------------------------------------------------------------------
bool MTDThresholdClusterizer::setup(const MTDGeometry* geom, const MTDTopology* topo, const DetId& id) 
{ 
  theCurrentId=id;
  //using geopraphicalId here
  const auto& thedet = geom->idToDet(id);
  if( thedet == nullptr ) {
    throw cms::Exception("MTDThresholdClusterizer") << "GeographicalID: " << std::hex
						    << id.rawId()
						    << " is invalid!" << std::dec
						    << std::endl;
    return false;
  }
  const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
  const RectangularMTDTopology& topol = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
  
  // Get the new sizes.
  unsigned int nrows = topol.nrows();      // rows in x
  unsigned int ncols = topol.ncolumns();   // cols in y
  
  theNumOfRows = nrows;  // Set new sizes
  theNumOfCols = ncols;
  
  DEBUG("Buffer size [" << theNumOfRows+1 << "," << theNumOfCols+1 << "]");
  
  if ( nrows > theBuffer.rows() || 
       ncols > theBuffer.columns() ) 
    { // change only when a larger is needed
      // Resize the buffer
      theBuffer.setSize(nrows+1,ncols+1);  // +1 needed for MTD
      bufferAlreadySet = true;
    }
  
  return true;   
}
//----------------------------------------------------------------------------
//!  \brief Cluster hits.
//!  This method operates on a matrix of hits
//!  and finds the largest contiguous cluster around
//!  each seed hit.
//!  Input and output data stored in DetSet
//----------------------------------------------------------------------------
void MTDThresholdClusterizer::clusterize( const FTLRecHitCollection & input,
					  const MTDGeometry* geom,
					  const MTDTopology* topo,
					  FTLClusterCollection& output) {
  
  FTLRecHitCollection::const_iterator begin = input.begin();
  FTLRecHitCollection::const_iterator end   = input.end();
  
  // Do not bother for empty detectors
  if (begin == end) 
    {
      edm::LogInfo("MTDThresholdClusterizer") << "No hits to clusterize";
      return;
    }

  DEBUG("Input collection " << input.size());
  assert(output.empty());

  std::set<unsigned> geoIds; 
  std::multimap<unsigned, unsigned> geoIdToIdx;
  
  unsigned index = 0;
  for(const auto& hit : input) 
    {
      MTDDetId mtdId=MTDDetId(hit.detid());
      if (mtdId.subDetector() != MTDDetId::FastTime)
	{
	  throw cms::Exception("MTDThresholdClusterizer") << "MTDDetId: " << std::hex
							  << mtdId.rawId()
							  << " is invalid!" << std::dec
							  << std::endl;
	}

      if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
	{
	  BTLDetId hitId(hit.detid());
	  DetId geoId = hitId.geographicalId( (BTLDetId::CrysLayout) topo->getMTDTopologyMode() ); //for BTL topology gives different layout id
	  geoIdToIdx.emplace(geoId,index);
	  geoIds.emplace(geoId);
	  ++index;
	}
      else if ( mtdId.mtdSubDetector() == MTDDetId::ETL )
	{
	  ETLDetId hitId(hit.detid());
	  DetId geoId = hitId.geographicalId();
	  geoIdToIdx.emplace(geoId,index);
	  geoIds.emplace(geoId);
	  ++index;
	}
      else
	{
	  throw cms::Exception("MTDThresholdClusterizer") << "MTDDetId: " << std::hex
							  << mtdId.rawId()
							  << " is invalid!" << std::dec
							  << std::endl;
	}
    }

  //cluster hits within geoIds (modules)
  for(unsigned id : geoIds) {
    //  Set up the clusterization on this DetId.
    if ( !setup(geom,topo,DetId(id)) ) 
      return;
    
    auto range = geoIdToIdx.equal_range(id);
    DEBUG("Matching Ids for " << std::hex << id << std::dec << " [" <<  range.first->second << "," << range.second->second << "]");
    
    //  Copy MTDRecHits to the buffer array; select the seed hits
    //  on the way, and store them in theSeeds.
    for(auto itr = range.first; itr != range.second; ++itr) {
      const unsigned hitidx = itr->second;
      copy_to_buffer(begin+hitidx, geom, topo);
    }
    
    FTLClusterCollection::FastFiller clustersOnDet(output,id);

    for (unsigned int i = 0; i < theSeeds.size(); i++) 
      {
	if ( theBuffer.energy(theSeeds[i]) > theSeedThreshold ) 
	  {  // Is this seed still valid?
	    //  Make a cluster around this seed
	    const FTLCluster & cluster = make_cluster( theSeeds[i] );
	    
	    //  Check if the cluster is above threshold  
	    if ( cluster.energy() > theClusterThreshold) 
	      {
		DEBUG("putting in this cluster " << i << " #hits:" << cluster.size() 
		      << " E:" << cluster.energy() 
		      << " T:" << cluster.time() 
		      << " X:" << cluster.x() 
		      << " Y:" << cluster.y());
		clustersOnDet.push_back( cluster ); 
	      }
	  }
      }
  
    // Erase the seeds.
    theSeeds.clear();  
    //  Need to clean unused hits from the buffer array.
    for(auto itr = range.first; itr != range.second; ++itr) {
      const unsigned hitidx = itr->second;
      clear_buffer(begin+hitidx);
    }
  }
}

//----------------------------------------------------------------------------
//!  \brief Clear the internal buffer array.
//!
//!  MTDs which are not part of recognized clusters are NOT ERASED 
//!  during the cluster finding.  Erase them now.
//!
//----------------------------------------------------------------------------
void MTDThresholdClusterizer::clear_buffer( RecHitIterator itr ) 
{
  theBuffer.clear( itr->row(), itr->column() ); 
}

//----------------------------------------------------------------------------
//! \brief Copy FTLRecHit into the buffer, identify seeds.
//----------------------------------------------------------------------------
void MTDThresholdClusterizer::copy_to_buffer( RecHitIterator itr, const MTDGeometry* geom, const MTDTopology* topo ) 
{
    MTDDetId mtdId=MTDDetId(itr->detid());
    int row = itr->row();
    int col = itr->column();
    int SubDet = 0;
    float energy = itr->energy();
    float time = itr->time();
    float timeError = itr->timeError();
    float position = itr->position().second;
    //float positionError = itr->positionError();
    //int new_col = round((position*64.)/57.);
    //float position = 0.0;
    //float positionError = 0.0;
    float x = -999.0; 
    float y = -999.0;
    float z = -999.0;
    float xError = 0.3;
    float yError = 0.3;		//error in cm
    float zError = 0.6;
    if( mtdId.mtdSubDetector() == MTDDetId::BTL )
    {
	  SubDet = 1;
          BTLDetId id = itr->id();
	  DetId geoId = id.geographicalId( (BTLDetId::CrysLayout) topo->getMTDTopologyMode() );
	  const auto& det = geom -> idToDet(geoId);
	  const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
	  const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
          MeasurementPoint mp(row,col);
	  LocalPoint lp_ctr = topo.localPosition(mp);
	  float z_offset = 2.85-(position/10.);  //need to add the offset as the distance from centre to edge (z direction in global and y direction in local)
	  LocalPoint lp(lp_ctr.x(),lp_ctr.y()+z_offset,lp_ctr.z());
	  //std::cout << "lp.x: " << lp.x() << " lp.y: " << lp.y() <<" pos: "<<position << std::endl;
	  GlobalPoint gp = det->toGlobal(lp);
  	  x = gp.x();
	  y = gp.y();
	  z = gp.z();
	  //float z_offset = 2.85-position;  //need to add the offset as the distance from centre to edge (z direction)
    }
    else if( mtdId.mtdSubDetector() == MTDDetId::ETL )
    {
	  SubDet = -1;
	  ETLDetId id = itr->id();
	  DetId geoId = id.geographicalId( );
	  const auto& det = geom -> idToDet(geoId);
	  const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
	  const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
	  
	  MeasurementPoint mp(row,col);
	  LocalPoint lp = topo.localPosition(mp);
	  GlobalPoint gp = det->toGlobal(lp);
	  x = gp.x();
	  y = gp.y();
    	  z = gp.z();
    }

    DEBUG("ROW " <<  row << " COL " << col << " ENERGY " << energy << " TIME " << time);
    if ( energy > theHitThreshold) {
      theBuffer.set( row, col, SubDet, energy, time, timeError, x, xError, y, yError, z, zError); 
      if ( energy > theSeedThreshold) theSeeds.push_back( FTLCluster::FTLHitPos(row,col));
      //sort seeds?
    }
}


//----------------------------------------------------------------------------
//!  \brief The actual clustering algorithm: group the neighboring hits around the seed.
//----------------------------------------------------------------------------
FTLCluster 
MTDThresholdClusterizer::make_cluster( const FTLCluster::FTLHitPos& hit ) 
{
  
  //First we acquire the seeds for the clusters
  const float thePositionThreshold = 0;   //in cm 

  float seed_subdet= theBuffer.SubDet(hit.row(),hit.col());
  float seed_energy= theBuffer.energy(hit.row(), hit.col());
  float seed_time= theBuffer.time(hit.row(), hit.col());
  float seed_time_error= theBuffer.time_error(hit.row(), hit.col());
  float seed_x= theBuffer.x(hit.row(), hit.col());
  float seed_x_error= theBuffer.x_error(hit.row(), hit.col());
  float seed_y= theBuffer.y(hit.row(), hit.col());
  float seed_y_error= theBuffer.y_error(hit.row(), hit.col());
  float seed_z= theBuffer.z(hit.row(), hit.col());
  float seed_z_error= theBuffer.z_error(hit.row(), hit.col());
  theBuffer.clear(hit);
  
  //std::cout << "seed x: " << seed_x <<"+-"<< seed_x_error << "  seed y: " << seed_y << "+-" << seed_y_error << std::endl;
  AccretionCluster acluster;
  acluster.add(hit, seed_energy, seed_time, seed_time_error);
  
  bool stopClus=false;
  //Here we search all hits adjacent to all hits in the cluster.
  while ( ! acluster.empty() && ! stopClus) 
    {
      //This is the standard algorithm to find and add a hit
      auto curInd = acluster.top(); acluster.pop();
      for ( auto c = std::max(0,int(acluster.y[curInd]-1)); c < std::min(int(acluster.y[curInd]+2),int(theBuffer.columns())) && !stopClus; ++c) {
	  for ( auto r = std::max(0,int(acluster.x[curInd]-1)); r < std::min(int(acluster.x[curInd]+2),int(theBuffer.rows()))  && !stopClus; ++r)  {
	  if ( theBuffer.energy(r,c) > theHitThreshold) {
	    if (std::abs(theBuffer.time(r,c) - seed_time) > theTimeThreshold*sqrt( theBuffer.time_error(r,c)*theBuffer.time_error(r,c) + seed_time_error*seed_time_error))
	      continue;
            if ((seed_subdet == 1) && (theBuffer.SubDet(r,c) == 1)){
	      if (sqrt(pow((theBuffer.x(r,c) - seed_x),2)+pow((theBuffer.y(r,c)-seed_y),2)+pow((theBuffer.z(r,c)-seed_z),2)) > thePositionThreshold*sqrt( theBuffer.x_error(r,c)*theBuffer.x_error(r,c) + seed_x_error*seed_x_error + theBuffer.y_error(r,c)*theBuffer.y_error(r,c) + seed_y_error*seed_y_error + theBuffer.z_error(r,c)*theBuffer.z_error(r,c) + seed_z_error*seed_z_error ))
                continue;
	    }
	    
	    FTLCluster::FTLHitPos newhit(r,c);
	    if (!acluster.add( newhit, theBuffer.energy(r,c), theBuffer.time(r,c), theBuffer.time_error(r,c))) 
	      {
		stopClus=true;
		break;
	      }
	    theBuffer.clear(newhit);
	  }
	}
      }
    }  // while accretion
  
  FTLCluster cluster( theCurrentId, acluster.isize, 
		      acluster.energy.data(), acluster.time.data(), acluster.timeError.data(), 
		      acluster.x.data(), acluster.y.data(), 
		      acluster.xmin, acluster.ymin);
  return cluster;
}
