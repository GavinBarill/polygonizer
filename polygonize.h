#ifndef IGL_POLYGONIZE_H
#define IGL_POLYGONIZE_H
#include <Eigen/Core>

namespace igl
{
  // Inputs:
  //   implicit  implicit function so that implicit(x,y,z)<0 means inside,
  //     implicit(x,y,z)>0 means outside and implicit(x,y,z)==0 means on the
  //     surface
  //   width  width of partitioning cubes
  //   range  maximum range of cubes (+/- on the three axes) from first cube
  //   seed  3d starting point on or near the surface
  // Outputs:
  //   V  #V by 3 list of vertex positions
  //   F  #F by 3 list of triangle indices
  //   N  #V by 3 list of vertex normals
  // Returns true iff exitted successfully
  template <
    typename Derivedseed, 
    typename DerivedV, 
    typename DerivedF,
    typename DerivedN>
  bool polygonize(
    const std::function< double(double,double,double) > & implicit,
    const double width,
    const int range,
    const Eigen::MatrixBase<Derivedseed> & seed,
    Eigen::PlainObjectBase<DerivedV> & V,
    Eigen::PlainObjectBase<DerivedF> & F,
    Eigen::PlainObjectBase<DerivedN> & N);
}

// Implementation

#include "Polygonizer.hpp"
#include <igl/list_to_matrix.h>
#include <vector>
template <
  typename Derivedseed, 
  typename DerivedV, 
  typename DerivedF,
  typename DerivedN>
bool igl::polygonize(
  const std::function< double(double,double,double) > & implicit,
  const double width,
  const int range,
  const Eigen::MatrixBase<Derivedseed> & seed,
  Eigen::PlainObjectBase<DerivedV> & V,
  Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedN> & N)
{
  std::vector<std::vector<typename DerivedV::Scalar > > vV,vN;
  std::vector<std::vector<typename DerivedF::Scalar > > vF;
  const auto & add_triangle = 
    [&vV,&vF,&vN]
    (int i1, int i2, int i3, std::vector<Polygonizer::VERTEX> & vertices)
      ->bool
  {
    // Add this triangle
    vF.emplace_back(std::vector<typename DerivedF::Scalar >{i3,i2,i1});
    // Add new vertices as they show up in list
    if(vertices.size() > vV.size())
    {
      for(int i = vV.size();i<vertices.size();i++)
      {
        const Polygonizer::VERTEX & v = vertices[i];
        vV.emplace_back(
          std::vector<typename DerivedV::Scalar >
          { v.position.x, v.position.y, v.position.z});
        vN.emplace_back(
          std::vector<typename DerivedV::Scalar >
          { v.normal.x, v.normal.y, v.normal.z});
      }
    }
    assert(vertices.size() == vV.size());
    return true;
  };
  std::string err = Polygonizer::polygonize(
    implicit,width,range,seed(0),seed(1),seed(2), 
    add_triangle,0);
  if(!err.empty())
  {
    return false;
  }
  igl::list_to_matrix(vV,V);
  igl::list_to_matrix(vF,F);
  igl::list_to_matrix(vN,N);
  return true;
}

#endif
