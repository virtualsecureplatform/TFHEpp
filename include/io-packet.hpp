#pragma once

#include "../thirdparties/cereal/include/cereal/archives/portable_binary.hpp"
#include "../thirdparties/cereal/include/cereal/cereal.hpp"
#include "../thirdparties/cereal/include/cereal/types/unordered_map.hpp"
#include "../thirdparties/cereal/include/cereal/types/string.hpp"
#include "../thirdparties/cereal/include/cereal/types/vector.hpp"
#include "../thirdparties/cereal/include/cereal/types/array.hpp"
#include "../thirdparties/cereal/include/cereal/types/optional.hpp"

#include "params.hpp"
#include <unordered_map>
#include <string>
#include <vector>

namespace TFHEpp {

	struct IOpacket{
		std::unordered_map<std::string,std::vector<TFHEpp::TLWE<TFHEpp::lvl0param>>> tlwelvl0ports;
		std::unordered_map<std::string,std::vector<TFHEpp::TRLWE<TFHEpp::lvl2param>>> trlwelvl2ports;
		std::unordered_map<std::string,std::vector<TFHEpp::TRGSW<TFHEpp::lvl2param>>> trgswlvl2ports;
		
		template <class Archive>
		void serialize(Archive &archive)
		{
			archive(tlwelvl0ports,trlwelvl2ports,trgswlvl2ports);
		}
	};
}