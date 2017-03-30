#include "../headers/ssp.h"
#include <algorithm>

HAC_USING_LIBRARY_NAMESPACE;

sspmod::sspmod()
{
	sonic_l.r = -1;
	sonic_l.h = -1;
	sonic_u.h = -1;
	sonic_u.r = -1;
	z_up = -1;
	z_low = -1;
	dist_number = dist_numbers::one;
};

ss_at_point sspmod::get_ssp(const float &_r, const float &_z)
{
	ss_at_point ss_at_point_out;
	memset(&ss_at_point_out, 0, sizeof(ss_at_point));

	if (dist_number == dist_numbers::one)
	{
		if (_z < z_low || _z > z_up)
		{
			Layer = depth_ssp.upper_bound(_z);
			z_up = Layer->second.h;
			Layer--;
			sonic_l = Layer->second;
			z_low = sonic_l.h;
		}

		ss_at_point_out.c = sonic_l.c + (_z - sonic_l.h) * sonic_l.cz;
//		ss_at_point_out.grad_c.r = 0.f;
		ss_at_point_out.grad_c.z = sonic_l.cz;

//		ss_at_point_out.crr = 0.f;
//		ss_at_point_out.crz = 0.f;
//		ss_at_point_out.czz = 0.f;
	}
	else
	{
		if (_r < sonic_l.r || _r > sonic_u.r)
		{
			iter_range_depth_ssp_up = range_depth_ssp.upper_bound(_r);
			iter_range_depth_ssp_low = iter_range_depth_ssp_up;
			--iter_range_depth_ssp_low;
			auto iter2 = iter_range_depth_ssp_up->second.upper_bound(_z);
			z_up = iter2->second.h;
			sonic_u = (--iter2)->second;
			iter2 = iter_range_depth_ssp_low->second.upper_bound(_z);
			sonic_l = (--iter2)->second;
			z_low = sonic_l.h;
		} 
		else
		{
			if (_z < z_low || _z > z_up)
			{
				auto iter2 = iter_range_depth_ssp_up->second.upper_bound(_z);

				z_up = iter2->second.h;
				sonic_u = (--iter2)->second;
				iter2 = iter_range_depth_ssp_low->second.upper_bound(_z);
				sonic_l = (--iter2)->second;
				z_low = sonic_l.h;
			}
		}

		// for this depth, _z get the sound speed at both ends of the segment
		auto c1 = sonic_l.c + (_z - sonic_l.h) * sonic_l.cz;
		auto c2 = sonic_u.c + (_z - sonic_u.h) * sonic_l.cz;

		// s = proportional distance of _r in range
		auto s = (_r - sonic_l.r) / (sonic_u.r - sonic_l.r);

		ss_at_point_out.c = (1.0f - s) * c1 + s * c2;
		ss_at_point_out.grad_c.r = (c2 - c1) / (sonic_u.r - sonic_l.r);
		ss_at_point_out.grad_c.z = (1.0f - s) * sonic_l.cz + s * sonic_u.cz;

		ss_at_point_out.crr = 0.f;
		ss_at_point_out.crz = 0.f;
		ss_at_point_out.czz = 0.f;
	}

	return ss_at_point_out;
};

error const sspmod::set_ssp(std::map<float, GDB::det::SSP> const &ssp_slices)
{
	error ret = error::OK_terminate;

	if (ssp_slices.empty())
	{
		ret = error::spp_param_not_set;
	}
	// ВРСЗ
	{
		sonic _sonic;
		if (ssp_slices.size() > 1)	// ВРСЗ меняется
		{
			dist_number = sspmod::dist_numbers::some;

			for (auto iter_r = ssp_slices.begin(); iter_r != ssp_slices.end(); iter_r++)
			{
				_sonic.r = iter_r->first;
				_sonic.c = iter_r->second.sonic[0];
				_sonic.h = iter_r->second.depth[0] - 1;
				std::map<float, sonic> depth_ssp;
				depth_ssp.insert({ _sonic.h, _sonic });
				for (size_t i_z = 0; i_z < iter_r->second.sonic.size(); i_z++)
				{
					_sonic.c = iter_r->second.sonic[i_z];
					_sonic.h = iter_r->second.depth[i_z];
					depth_ssp.insert({ _sonic.h, _sonic });
				}
				_sonic.h += 100;
				depth_ssp.insert({ _sonic.h, _sonic });

				for (auto z_iter = depth_ssp.begin(); z_iter != --depth_ssp.end(); z_iter++)
				{
					auto z_iter_next = z_iter; z_iter_next++;
					z_iter->second.cz = (z_iter_next->second.c - z_iter->second.c) / (z_iter_next->second.h - z_iter->second.h);
				}

				range_depth_ssp.insert({ _sonic.r, depth_ssp });
			}
			// учет возможых выходов за пределы расчетов
			{
				auto iter_ssp_high = --range_depth_ssp.end();
				auto ssp_high = iter_ssp_high->second;
				auto r = iter_ssp_high->first * 1.1f;
				for (auto iter_s = ssp_high.begin(); iter_s != ssp_high.end(); iter_s++)
				{
					iter_s->second.r = r;
				}
				range_depth_ssp.insert({ r, ssp_high });

				auto ssp_low = range_depth_ssp.begin()->second;
				for (auto iter_s = ssp_low.begin(); iter_s != ssp_low.end(); iter_s++)
				{
					iter_s->second.r = -r;
				}
				range_depth_ssp.insert({ -r, ssp_low });
			}
		}
		else // ВРСЗ неизменно
		{
			dist_number = sspmod::dist_numbers::one;
			for (auto iter_r = ssp_slices.begin(); iter_r != ssp_slices.end(); iter_r++)
			{
				_sonic.r = iter_r->first;
				_sonic.c = iter_r->second.sonic[0];
				_sonic.h = iter_r->second.depth[0] - 1;
				depth_ssp.insert({ _sonic.h, _sonic });
				for (size_t i_z = 0; i_z < iter_r->second.sonic.size(); i_z++)
				{
					_sonic.c = iter_r->second.sonic[i_z];
					_sonic.h = iter_r->second.depth[i_z];
					depth_ssp.insert({ _sonic.h, _sonic });
				}
				_sonic.h += 100;
				depth_ssp.insert({ _sonic.h, _sonic });

				for (auto z_iter = depth_ssp.begin(); z_iter != --depth_ssp.end(); z_iter++)
				{
					auto z_iter_next = z_iter; z_iter_next++;
					z_iter->second.cz = (z_iter_next->second.c - z_iter->second.c) / (z_iter_next->second.h - z_iter->second.h);
				}
			}
		}
	}
	
	return ret;
};
