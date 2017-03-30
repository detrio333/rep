#include "../headers/field_reverb.h"

HAC_USING_LIBRARY_NAMESPACE;

const float field_reverb::get_pwr(time_band const&_time, space_band_v const&_space) const
{
	float P_rev = 0;

	for (auto ray = Q_angle_time.upper_bound(_space.low); ray != --Q_angle_time.lower_bound(_space.high); ray++)
	{
		P_rev += ray->second.get_pwr(_time);
	}

// 	for (auto ray : Q_angle_time)
// 	{
// 		if (_space.angle_contains(ray.first))
// 		{
// 			P_rev += ray.second.get_pwr(_time);
// 		}
// 	}

	return P_rev;
};

const float field_reverb::field_reverb_slice::get_pwr(time_band const&_time) const
{
	float P_rev = 0;

	for (auto time = Q_time.upper_bound(_time.begin); time != --Q_time.lower_bound(_time.end); time++)
	{
		P_rev += time->second;
	}

// 	for (auto time : Q_time)
// 	{
// 		auto t = time.first;
// 		if (t > _time.begin && t <_time.end)
// 		{
// 			P_rev += time.second;
// 		}
// 	}

	return P_rev;
};