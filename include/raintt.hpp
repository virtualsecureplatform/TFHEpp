#pragma once

namespace raintt{
	template <typename T>
	constexpr T ipow(T num, unsigned int pow)
	{
		return (pow >= sizeof(unsigned int)*8) ? 0 :
			pow == 0 ? 1 : num * ipow(num, pow-1);
	}
	constexpr uint32_t k = 5;
	constexpr uint32_t P = ((ipow<uint32_t>(5,4))ULL<<20)+1;
	constexpr uint32_t W = 11;

	
}	//namespace raintt