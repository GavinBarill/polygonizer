#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
export PS4=""
set -o xtrace
source ~/.profile
echo $DYLD_FALLBACK_LIBRARY_PATH
TEMP="$0.cpp"
printf "//" | cat - $0 >$TEMP
#/usr/local/opt/llvm/bin/clang++ -fopenmp=libiomp5 -I /usr/local/include/libiomp/ -liomp5 -g -O0 -std=c++11 -ferror-limit=4 -o .main $TEMP -msse4.2 \
#g++-4.8   -g -O0 -std=c++11 -fopenmp -o .main $TEMP -DNDEBUG -msse4.2 -m32 \
export LIBIGL=/usr/local/libigl/
#export LIBIGL=/Users/ajx/Dropbox/interactive-segmentation/libigl/
#-DIGL_STATIC_LIBRARY -L"$LIBIGL"/lib -ligl_viewer -ligl -ligl_matlab -ligl_opengl -ligl_opengl_glfw -ligl_cgal -L"$LIBIGL"/lib -lglfw3 \
#clang++ -g -O0 -std=c++11 -o .main $TEMP -msse4.2 \
clang++ -g -O3 -std=c++11 -o .main $TEMP -DNDEBUG -msse4.2 \
  -I. \
  -L"$LIBIGL"/lib  -lglfw3 \
  -I"$LIBIGL"/external/nanogui/ext/eigen/ \
  -I"$LIBIGL"/include \
  -I"$LIBIGL"/external/AntTweakBar/include \
  -I"$LIBIGL"/external/AntTweakBar/src \
  -I"$LIBIGL"/external/tetgen \
  -I"$LIBIGL"/external/tinyxml2/ \
  -I"$LIBIGL"/external/Singular_Value_Decomposition/ \
  -framework Carbon -framework QuartzCore -framework IOKit \
  -I"$LIBIGL"/external/nanogui/include \
  -I"$LIBIGL"/external/nanogui/ext \
  -I"$LIBIGL"/external/nanogui/ext/nanovg/src \
  -I"$LIBIGL"/external/embree/ \
  -I"$LIBIGL"/external/embree/include \
  -L/usr/local/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt \
  -framework OpenGL \
  -framework AppKit \
  -I/Applications/MATLAB_R2017a.app/extern/include/ \
  -L/Applications/MATLAB_R2017a.app/bin/maci64/ -lmx -lmat -lmex -leng \
  -L/usr/local/lib -lboost_thread-mt -lboost_system-mt \
  -L/usr/local/lib -lboost_program_options-mt -fno-math-errno \
  -I/Users/ajx/Dropbox/ \
  -I/usr/local/libigl/external/nanogui/ext/glfw/include \
&& /Applications/Xcode.app/Contents/Developer/usr/bin/lldb -b -o r ./.main -- "$@"
#&& ./.main "$@"
#-I~/Documents/eigen/ \
#-framework GLUT \
#  -L/usr/local/libigl/external/glfw/lib -lglfw3 -framework Carbon -framework QuartzCore -framework IOKit \
#  -I/usr/local/libigl/external/nanogui/include \
#  -I/usr/local/libigl/external/nanogui/ext \
#  -I/usr/local/libigl/external/nanogui/ext/nanovg/src \
#  -L/usr/local/libigl/external/nanogui/build -lnanogui \
#  -L/usr/local/libigl/external/nanogui/ext/nanovg \
#  -L/usr/local/libigl/lib/libigl/embree -lembree \
#  -I/usr/local/mosek/7/tools/platform/osx64x86/h \
#  -L/usr/local/mosek/7/tools/platform/osx64x86/bin -lmosek64 \
#rm -f .main
#rm -f $TEMP
# END BASH SCRIPT
exit
*/

#include "polygonize.h"
#include <igl/writeOBJ.h>
#include <Eigen/Core>
int main(int argc, char * argv[])
{
  const auto & sphere = 
    [](const double x, const double y ,const double z)->double
  {
    return sqrt(x*x+y*y+z*z)-1.0;
  };
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  igl::polygonize(sphere,0.1,20,Eigen::RowVector3d(0,0,0),V,F,N);
  igl::writeOBJ("test.obj",V,F,N,F,V,F);
}

