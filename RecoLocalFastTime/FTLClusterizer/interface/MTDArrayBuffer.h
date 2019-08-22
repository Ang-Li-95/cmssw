#ifndef RecoLocalTracker_FTLClusterizer_MTDArrayBuffer_H
#define RecoLocalTracker_FTLClusterizer_MTDArrayBuffer_H

//----------------------------------------------------------------------------
//! \class MTDArrayBuffer
//! \brief Class to store ADC counts and times during clustering.
//!
//----------------------------------------------------------------------------

// We use FTLHitPos which is an inner class of FTLCluster:
#include "DataFormats/FTLRecHit/interface/FTLCluster.h"

#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include <vector>
#include <iostream>

class MTDArrayBuffer 
{
 public:
  typedef unsigned int uint;
  
  inline MTDArrayBuffer( uint rows, uint cols);
  inline MTDArrayBuffer( ){}
  
  inline void setSize( uint rows, uint cols);

  //add SubDet ti identify whether the Hit is in BTL(+1) or ETL(-1), if not clear then use 0
  inline int SubDet( uint row, uint col) const;
  inline int SubDet( const FTLCluster::FTLHitPos&) const;
  inline float energy( uint row, uint col) const;
  inline float energy( const FTLCluster::FTLHitPos&) const;
  inline float time( uint row, uint col) const;
  inline float time( const FTLCluster::FTLHitPos&) const;
  inline float time_error( uint row, uint col) const;
  inline float time_error( const FTLCluster::FTLHitPos&) const;
  
  inline LocalError local_error( uint row, uint col) const;
  inline LocalError local_error( const FTLCluster::FTLHitPos&) const;
  inline GlobalPoint global_point( uint row, uint col) const;
  inline GlobalPoint global_point( const FTLCluster::FTLHitPos&) const;
  
  //inline float x( uint row, uint col) const;
  //inline float x( const FTLCluster::FTLHitPos&) const;
  //inline float x_error( uint row, uint col) const;
  //inline float x_error( const FTLCluster::FTLHitPos&) const;
  //inline float y( uint row, uint col) const;
  //inline float y( const FTLCluster::FTLHitPos&) const;
  //inline float y_error( uint row, uint col) const;
  //inline float y_error( const FTLCluster::FTLHitPos&) const;
  //inline float z( uint row, uint col) const;
  //inline float z( const FTLCluster::FTLHitPos&) const;
  //inline float z_error( uint row, uint col) const;
  //inline float z_error( const FTLCluster::FTLHitPos&) const;


  inline uint rows() const { return nrows;}
  inline uint columns() const { return ncols;}

  inline bool inside(uint row, uint col) const;

  inline void clear(uint row, uint col) 
  {
    LocalError le_n(0,0,0);
    GlobalPoint gp_n(0,0,0);
    set_SubDet( row, col, 0);
    set_energy( row, col, 0.);
    set_time( row, col, 0.);
    set_time_error( row, col, 0.);
    set_local_error( row, col, le_n);
    set_global_point( row, col, gp_n);
    //set_x( row, col, 0.);
    //set_x_error( row, col, 0.);
    //set_y( row, col, 0.);
    //set_y_error( row, col, 0.);
    //set_z( row, col, 0.);
    //set_z_error( row, col, 0.);

  }
  inline void clear(const FTLCluster::FTLHitPos& pos) 
  {
    clear(pos.row(),pos.col());
  }

  inline void set( uint row, uint col, float energy, float time, float time_error);
  inline void set( const FTLCluster::FTLHitPos&, float energy, float time, float time_error);
  inline void set( uint row, uint col, int SubDet, float energy, float time, float time_error, const LocalError& local_error, const GlobalPoint& global_point);
  inline void set( const FTLCluster::FTLHitPos&, int SubDet, float energy, float time, float time_error, const LocalError& local_error, const GlobalPoint& global_point);


  inline void set_SubDet( uint row, uint col, int SubDet);
  inline void set_SubDet( const FTLCluster::FTLHitPos&, int SubDet);

  inline void set_energy( uint row, uint col, float energy);
  inline void set_energy( const FTLCluster::FTLHitPos&, float energy);
  inline void add_energy( uint row, uint col, float energy);

  inline void set_time( uint row, uint col, float time);
  inline void set_time( const FTLCluster::FTLHitPos&, float time);

  inline void set_time_error( uint row, uint col, float time_error);
  inline void set_time_error( const FTLCluster::FTLHitPos&, float time_error);

  inline void set_global_point( uint row, uint col, const GlobalPoint& gp);
  inline void set_global_point( const FTLCluster::FTLHitPos&, const GlobalPoint& gp); 

  inline void set_local_error( uint row, uint col, const LocalError& le);
  inline void set_local_error( const FTLCluster::FTLHitPos&, const LocalError& le);

  //inline void set_x( uint row, uint col, float x);
  //inline void set_x( const FTLCluster::FTLHitPos&, float x);

  //inline void set_x_error( uint row, uint col, float x_error);
  //inline void set_x_error( const FTLCluster::FTLHitPos&, float x_error);

  //inline void set_y( uint row, uint col, float y);
  //inline void set_y( const FTLCluster::FTLHitPos&, float y);

  //inline void set_y_error( uint row, uint col, float y_error);
  //inline void set_y_error( const FTLCluster::FTLHitPos&, float y_error);

  //inline void set_z( uint row, uint col, float z);
  //inline void set_z( const FTLCluster::FTLHitPos&, float z);

  //inline void set_z_error( uint row, uint col, float z_error);
  //inline void set_z_error( const FTLCluster::FTLHitPos&, float z_error);


  uint size() const { return hitEnergy_vec.size();}

  /// Definition of indexing within the buffer.
  uint index( uint row, uint col) const {return col*nrows+row;}
  uint index( const FTLCluster::FTLHitPos& pix) const { return index(pix.row(), pix.col()); }

 private:
  std::vector<int> hitSubDet_vec;
  std::vector<float> hitEnergy_vec;   
  std::vector<float> hitTime_vec;  
  std::vector<float> hitTimeError_vec;   
  std::vector<GlobalPoint> hitGP_vec;
  std::vector<LocalError> hitLE_vec;
  //std::vector<float> hitX_vec;
  //std::vector<float> hitXError_vec;
  //std::vector<float> hitY_vec;
  //std::vector<float> hitYError_vec;
  //std::vector<float> hitZ_vec;
  //std::vector<float> hitZError_vec;
  uint nrows;
  uint ncols;
};

MTDArrayBuffer::MTDArrayBuffer( uint rows, uint cols) 
  : hitSubDet_vec(rows*cols,0), hitEnergy_vec(rows*cols,0), hitTime_vec(rows*cols,0), hitTimeError_vec(rows*cols,0), hitGP_vec(rows*cols), hitLE_vec(rows*cols), nrows(rows), ncols(cols) {}

void MTDArrayBuffer::setSize( uint rows, uint cols) {
  hitSubDet_vec.resize(rows*cols,0);
  hitEnergy_vec.resize(rows*cols,0);
  hitTime_vec.resize(rows*cols,0);
  hitTimeError_vec.resize(rows*cols,0);
  hitGP_vec.resize(rows*cols);
  hitLE_vec.resize(rows*cols);
  //hitX_vec.resize(rows*cols,0);
  //hitXError_vec.resize(rows*cols,0);
  //hitY_vec.resize(rows*cols,0);
  //hitYError_vec.resize(rows*cols,0);
  //hitZ_vec.resize(rows*cols,0);
  //hitZError_vec.resize(rows*cols,0);
  nrows = rows;
  ncols = cols;
}

bool MTDArrayBuffer::inside(uint row, uint col) const 
{
  return ( row < nrows && col < ncols);
}

int MTDArrayBuffer::SubDet(uint row, uint col) const { return hitSubDet_vec[index(row,col)];}
int MTDArrayBuffer::SubDet(const FTLCluster::FTLHitPos& pix) const {return hitSubDet_vec[index(pix)];}

float MTDArrayBuffer::energy(uint row, uint col) const { return hitEnergy_vec[index(row,col)];}
float MTDArrayBuffer::energy(const FTLCluster::FTLHitPos& pix) const {return hitEnergy_vec[index(pix)];}

float MTDArrayBuffer::time(uint row, uint col) const  { return hitTime_vec[index(row,col)];}
float MTDArrayBuffer::time(const FTLCluster::FTLHitPos& pix) const {return hitTime_vec[index(pix)];}

float MTDArrayBuffer::time_error(uint row, uint col) const  { return hitTimeError_vec[index(row,col)];}
float MTDArrayBuffer::time_error(const FTLCluster::FTLHitPos& pix) const {return hitTimeError_vec[index(pix)];}

LocalError MTDArrayBuffer::local_error(uint row, uint col) const {return hitLE_vec[index(row,col)];}
LocalError MTDArrayBuffer::local_error(const FTLCluster::FTLHitPos& pix) const {return hitLE_vec[index(pix)];}

GlobalPoint MTDArrayBuffer::global_point(uint row, uint col) const {return hitGP_vec[index(row,col)];}
GlobalPoint MTDArrayBuffer::global_point(const FTLCluster::FTLHitPos& pix) const {return hitGP_vec[index(pix)];}

//float MTDArrayBuffer::x(uint row, uint col) const  { return hitX_vec[index(row,col)];}
//float MTDArrayBuffer::x(const FTLCluster::FTLHitPos& pix) const {return hitX_vec[index(pix)];}
//
//float MTDArrayBuffer::x_error(uint row, uint col) const  { return hitXError_vec[index(row,col)];}
//float MTDArrayBuffer::x_error(const FTLCluster::FTLHitPos& pix) const {return hitXError_vec[index(pix)];}
//
//float MTDArrayBuffer::y(uint row, uint col) const  { return hitY_vec[index(row,col)];}
//float MTDArrayBuffer::y(const FTLCluster::FTLHitPos& pix) const {return hitY_vec[index(pix)];}
//
//float MTDArrayBuffer::y_error(uint row, uint col) const  { return hitYError_vec[index(row,col)];}
//float MTDArrayBuffer::y_error(const FTLCluster::FTLHitPos& pix) const {return hitYError_vec[index(pix)];}
//
//float MTDArrayBuffer::z(uint row, uint col) const  { return hitZ_vec[index(row,col)];}
//float MTDArrayBuffer::z(const FTLCluster::FTLHitPos& pix) const {return hitZ_vec[index(pix)];}
//
//float MTDArrayBuffer::z_error(uint row, uint col) const  { return hitZError_vec[index(row,col)];}
//float MTDArrayBuffer::z_error(const FTLCluster::FTLHitPos& pix) const {return hitZError_vec[index(pix)];}

void MTDArrayBuffer::set( uint row, uint col, float energy, float time, float time_error) 
{
  hitEnergy_vec[index(row,col)] = energy;
  hitTime_vec[index(row,col)] = time;
  hitTimeError_vec[index(row,col)] = time_error;
}
void MTDArrayBuffer::set( const FTLCluster::FTLHitPos& pix, float energy, float time, float time_error) 
{
  set( pix.row(), pix.col(), energy, time, time_error);
}

void MTDArrayBuffer::set( uint row, uint col, int SubDet, float energy, float time, float time_error, const LocalError& local_error, const GlobalPoint& global_point)
{
  hitSubDet_vec[index(row,col)] = SubDet;
  hitEnergy_vec[index(row,col)] = energy;
  hitTime_vec[index(row,col)] = time;
  hitTimeError_vec[index(row,col)] = time_error;
  hitGP_vec[index(row,col)] = global_point;
  hitLE_vec[index(row,col)] = local_error;
  //hitX_vec[index(row,col)] = x;
  //hitXError_vec[index(row,col)] = x_error;
  //hitY_vec[index(row,col)] = y;
  //hitYError_vec[index(row,col)] = y_error;
  //hitZ_vec[index(row,col)] = z;
  //hitZError_vec[index(row,col)] = z_error;
  //std::cout << "vec_gp: " << hitGP_vec[index(row,col)].x() << " " << hitGP_vec[index(row,col)].y() << " " << hitGP_vec[index(row,col)].z() << std::endl;
  //std::cout << "vec_le: " << hitLE_vec[index(row,col)].xx() << " " << hitLE_vec[index(row,col)].yy() << std::endl;
}
void MTDArrayBuffer::set( const FTLCluster::FTLHitPos& pix, int SubDet, float energy, float time, float time_error, const LocalError& local_error, const GlobalPoint& global_point)
{
  set( pix.row(), pix.col(), SubDet, energy, time, time_error, local_error, global_point);
}

void MTDArrayBuffer::set_SubDet( uint row, uint col, int SubDet)
{
  hitSubDet_vec[index(row,col)] = SubDet;
}
void MTDArrayBuffer::set_SubDet( const FTLCluster::FTLHitPos& pix, int SubDet)
{
  hitSubDet_vec[index(pix)] = SubDet;
}

void MTDArrayBuffer::set_energy( uint row, uint col, float energy) 
{
  hitEnergy_vec[index(row,col)] = energy;
}
void MTDArrayBuffer::set_energy( const FTLCluster::FTLHitPos& pix, float energy)
{
  hitEnergy_vec[index(pix)] = energy;
}
void MTDArrayBuffer::add_energy( uint row, uint col, float energy)
{
  hitEnergy_vec[index(row,col)] += energy;
}

void MTDArrayBuffer::set_time( uint row, uint col, float time) 
{
  hitTime_vec[index(row,col)] = time;
}
void MTDArrayBuffer::set_time( const FTLCluster::FTLHitPos& pix, float time)
{
  hitTime_vec[index(pix)] = time;
}

void MTDArrayBuffer::set_time_error( uint row, uint col, float time_error) 
{
  hitTimeError_vec[index(row,col)] = time_error;
}
void MTDArrayBuffer::set_time_error( const FTLCluster::FTLHitPos& pix, float time_error)
{
  hitTimeError_vec[index(pix)] = time_error;
}

void MTDArrayBuffer::set_global_point( uint row, uint col, const GlobalPoint& gp)
{
  hitGP_vec[index(row,col)] = gp;
}
void MTDArrayBuffer::set_global_point( const FTLCluster::FTLHitPos& pix, const GlobalPoint& gp)
{
  hitGP_vec[index(pix)] = gp;
}

void MTDArrayBuffer::set_local_error( uint row, uint col, const LocalError& le)
{
  hitLE_vec[index(row,col)] = le;
}
void MTDArrayBuffer::set_local_error( const FTLCluster::FTLHitPos& pix, const LocalError& le)
{
  hitLE_vec[index(pix)] = le;
}

//void MTDArrayBuffer::set_x( uint row, uint col, float x)
//{
//  hitX_vec[index(row,col)] = x;
//}
//void MTDArrayBuffer::set_x( const FTLCluster::FTLHitPos& pix, float x)
//{
//  hitX_vec[index(pix)] = x;
//}
//
//void MTDArrayBuffer::set_x_error( uint row, uint col, float x_error)
//{
//  hitXError_vec[index(row,col)] = x_error;
//}
//void MTDArrayBuffer::set_x_error( const FTLCluster::FTLHitPos& pix, float x_error)
//{
//  hitXError_vec[index(pix)] = x_error;
//}
//
//void MTDArrayBuffer::set_y( uint row, uint col, float y)
//{
//  hitY_vec[index(row,col)] = y;
//}
//void MTDArrayBuffer::set_y( const FTLCluster::FTLHitPos& pix, float y)
//{
//  hitY_vec[index(pix)] = y;
//}
//
//void MTDArrayBuffer::set_y_error( uint row, uint col, float y_error)
//{
//  hitYError_vec[index(row,col)] = y_error;
//}
//void MTDArrayBuffer::set_y_error( const FTLCluster::FTLHitPos& pix, float y_error)
//{
//  hitYError_vec[index(pix)] = y_error;
//}
//
//void MTDArrayBuffer::set_z( uint row, uint col, float z)
//{
//  hitZ_vec[index(row,col)] = z;
//}
//void MTDArrayBuffer::set_z( const FTLCluster::FTLHitPos& pix, float z)
//{
//  hitZ_vec[index(pix)] = z;
//}
//
//void MTDArrayBuffer::set_z_error( uint row, uint col, float z_error)
//{
//  hitZError_vec[index(row,col)] = z_error;
//}
//void MTDArrayBuffer::set_z_error( const FTLCluster::FTLHitPos& pix, float z_error)
//{
//  hitZError_vec[index(pix)] = z_error;
//}
#endif
