#include "../headers/field_clutter.h"

HAC_USING_LIBRARY_NAMESPACE;

float const field_clutter::get_pwr(space_band_v const&_space_band_v)const
{
	float p = 0;

	for (auto angle_iter = Q.begin(); angle_iter != Q.end(); angle_iter++)
	{
		if (angle_iter->first >= _space_band_v.low && angle_iter->first <= _space_band_v.high)
		{
			p += angle_iter->second;
		}
	}

	return p*P;
};
