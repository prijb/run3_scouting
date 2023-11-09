// -*- C++ -*-
//Add includes for your classes here
#include <vector>
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
   struct dictionary {
      std::vector<std::vector<float> > vf2d;
      edm::Wrapper<std::vector<std::vector<float> > > wvf2d;

      std::vector<std::vector<bool> > vb2d;
      edm::Wrapper<std::vector<std::vector<bool> > > wvb2d;

      std::vector<std::vector<int> > vi2d;
      edm::Wrapper<std::vector<std::vector<int> > > wvi2d;

      std::vector<std::vector<std::vector<float> > > vf3d;
      edm::Wrapper<std::vector<std::vector<std::vector<float> > > > wvf3d;

      std::vector<std::vector<std::vector<bool> > > vb3d;
      edm::Wrapper<std::vector<std::vector<std::vector<bool> > > > wvb3d;

      std::vector<std::vector<std::vector<int> > > vi3d;
      edm::Wrapper<std::vector<std::vector<std::vector<int> > > > wvi3d;
      
   };
}
