#include "../headers/bands.h"

HAC_USING_LIBRARY_NAMESPACE;

freq_band::freq_band()
:low(std::numeric_limits<float>::quiet_NaN()), high(std::numeric_limits<float>::quiet_NaN()){

};

freq_band::freq_band(float const __low, float const __high)
:low(__low), high(__high){

};

space_band_v::space_band_v()
:low(std::numeric_limits<float>::quiet_NaN()), high(std::numeric_limits<float>::quiet_NaN()){

};

space_band_v::space_band_v(float const __low, float const __high)
:low(__low), high(__high){

};

time_band::time_band()
:begin(std::numeric_limits<int>::quiet_NaN()), end(std::numeric_limits<int>::quiet_NaN()){

}

time_band::time_band(unsigned int const __begin, unsigned int const __end)
: begin(__begin), end(__end){

};

float const freq_band::get_freq_avg()const
{
	return sqrt(low*high);
}

float const space_band_v::get_space_avg()const
{
	return (low + high) / 2;
};

const bool space_band_v::operator == (const space_band_v &_sb) const
{
	return (high == _sb.high && low == _sb.low);
};

const bool space_band_v::angle_contains(const float &_angle) const
{
	return (_angle <= high && _angle > low);
};