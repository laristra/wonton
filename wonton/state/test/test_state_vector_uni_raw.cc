/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/support/wonton.h"
#include "wonton/state/state_vector_uni_raw.h"

#include "gtest/gtest.h"

using namespace Wonton;

TEST(StateUni, BasicInt) {

	std::string name{"field"};
	
	StateVectorUniRaw<int> v(name);
  
  ASSERT_EQ(v.get_name(), name);
  ASSERT_EQ(v.get_type(), Field_type::MESH_FIELD);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(StateUni, BasicDouble1) {

	std::string name{"field"};
	
	StateVectorUniRaw<> v(name);
  
  ASSERT_EQ(v.get_name(), name);
  ASSERT_EQ(v.get_type(), Field_type::MESH_FIELD);
  ASSERT_EQ(v.data_type(), typeid(double));
  ASSERT_EQ(v.get_data(), nullptr);
  
}

TEST(StateUni, DataAccess) {

	std::string name{"field"};
	std::vector<double> data {1.,2.,3.};	
	
	StateVectorUniRaw<double> sv(name, Entity_kind::CELL, data.data());
	
	const double * sv_data{sv.get_data()};

	for (unsigned i=0; i < data.size(); i++) {
		ASSERT_EQ(sv_data[i],data[i]);
	}
}


