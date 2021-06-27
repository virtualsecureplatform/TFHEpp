#pragma once

#include "params.hpp"
#include <map>
#include <string>
#include <vector>

namespace TFHEpp {
	struct IOpacket{
		std::map<std::string,std::vector<TFHEpp::TLWE<TFHEpp::lvl0param>>> tlwelvl0ports;
		std::map<std::string,std::vector<TFHEpp::TRLWE<TFHEpp::lvl2param>>> trlwelvl2ports;
		std::map<std::string,std::vector<TFHEpp::TRGSW<TFHEpp::lvl2param>>> trgswlvl2ports;
		
		template <class Archive>
		void serialize(Archive &archive)
		{
			archive(tlwelvl0ports, trlwelvl2ports,trgswlvl2ports);
		}
	};
}