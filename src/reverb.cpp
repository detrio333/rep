#include "../headers/reverb.h"
#include <algorithm>

HAC_USING_LIBRARY_NAMESPACE;

error const reverb::set_param(calc_param const &_local_calc_param)
{
	error ret = error::OK_terminate;

	local_calc_param = _local_calc_param;

	auto result = find(local_calc_param.type.begin(), local_calc_param.type.end(), type_reverb::surface);
	result == local_calc_param.type.end() ? surface_r = false : surface_r = true;

	result = find(local_calc_param.type.begin(), local_calc_param.type.end(), type_reverb::volume);
	result == local_calc_param.type.end() ? volume_r = false : volume_r = true;

	result = find(local_calc_param.type.begin(), local_calc_param.type.end(), type_reverb::seabed);
	result == local_calc_param.type.end() ? seabed_r = false : seabed_r = true;

	return ret;
};

error const reverb::get_clutter(
	field_reverb &_field_reverb, const std::map<float, std::vector<step_ray>> &_rays_field
	)const
{
	error ret = error::OK_terminate;

	// если не пуст контейнер, то очищаем
	if (!_field_reverb.Q_angle_time.empty()) _field_reverb.Q_angle_time.clear();

	// дан контейнер с первым индексом по углу, и вторым по времени
	// надо выдать с первым по времени и вторым по углу (уже не надо)

	auto koef_attenuation = local_hydrology->koef_attenuation / 8685.8896f;	// в Nepers/m из dB/km
	auto dalpha = abs(_rays_field.begin()->first - (++_rays_field.begin())->first);

	auto src_coord = _rays_field.begin()->second.begin()->x;
//	local_hydrology->koef_attenuation = 0;

	for (auto ray : _rays_field)
	{
		field_reverb::field_reverb_slice _field_reverb_slice;

		bool end = false;

		for (auto _step_ray = ray.second.begin() + 1; _step_ray != ray.second.end(); _step_ray++)
		{
			if (_step_ray->tau * 2 > local_calc_param.t_lim)
			{
				end = true;
				break;
			}
			switch (_step_ray->event)
			{
			case step_ray::events::volume:
				if (volume_r)
				{
					float Wnorm_rev = koef_attenuation * _step_ray->Len * 2;
					auto x = src_coord.r - _step_ray->x.r;
					auto z = src_coord.z - _step_ray->x.z;
					auto dist = x*x + z*z;	// квадрат расстояния

					auto k_scatt = (--(--local_hydrology->scat.upper_bound(_step_ray->x.r))->second.upper_bound(_step_ray->x.z))->second;
					Wnorm_rev = local_calc_param.t_sign * 0.5f *_step_ray->c	* dalpha// рассеевающий объем
						* exp(-Wnorm_rev) / dist	// затухание при распространении
						* k_scatt	// коэф обратного рассеяния
						* _step_ray->Rfa*_step_ray->Rfa;	// потери при распространении	

					_field_reverb_slice.Q_time.insert({ _step_ray->tau * 2, Wnorm_rev });
				}
				break;
			case step_ray::events::surface:
				if (surface_r)
				{
					// угол скольжения
					auto angle = -atan2(_step_ray->Tray.z, _step_ray->Tray.r);
					
					float Wnorm_rev = koef_attenuation * _step_ray->Len * 2;
					auto x = src_coord.r - _step_ray->x.r;
					auto z = src_coord.z - _step_ray->x.z;
					auto dist = sqrt(x*x + z*z);	// квадрат расстояния

					// identify the top segment above the source
					auto top = (--local_hydrology->top_brdy.map_brdy.upper_bound(_step_ray->x.r))->second;

					auto sca_u = top.ScaSet.upper_bound(angle);
					auto sca_l = sca_u;
					sca_l--;

					float Alpha = (angle - sca_l->first) / (sca_u->first - sca_l->first);
					float k_scatt = (1 - Alpha) * sca_l->second + Alpha * sca_u->second;

					Wnorm_rev = local_calc_param.t_sign * 0.5f *_step_ray->c // рассеевающий объем
						* exp(-Wnorm_rev) / (dist*dist*dist)	// затухание при распространении
						* k_scatt	// коэф обратного рассеяния
						* _step_ray->Rfa*_step_ray->Rfa;	// потери при распространении

					_field_reverb_slice.Q_time.insert({ _step_ray->tau * 2, Wnorm_rev });
				}
				break;
			case step_ray::events::seabed:
				if (seabed_r)
				{
					// угол скольжения
					auto angle = -atan2(_step_ray->Tray.z, _step_ray->Tray.r);

					float Wnorm_rev = koef_attenuation * _step_ray->Len * 2;
					auto x = src_coord.r - _step_ray->x.r;
					auto z = src_coord.z - _step_ray->x.z;
					auto dist = sqrt(x*x + z*z);	// квадрат расстояния

					// identify the bottom segment above the source
					auto bot = (--local_hydrology->bot_brdy.map_brdy.upper_bound(_step_ray->x.r))->second;

					auto sca_u = bot.ScaSet.upper_bound(angle);
					auto sca_l = sca_u;
					sca_l--;

					float Alpha = (angle - sca_l->first) / (sca_u->first - sca_l->first);
					float k_scatt = (1 - Alpha) * sca_l->second + Alpha * sca_u->second;

					Wnorm_rev = local_calc_param.t_sign * 0.5f *_step_ray->c // рассеевающий объем
						* exp(-Wnorm_rev) / (dist*dist*dist)	// затухание при распространении
						* k_scatt	// коэф обратного рассеяния
						* _step_ray->Rfa*_step_ray->Rfa;	// потери при распространении

					_field_reverb_slice.Q_time.insert({ _step_ray->tau * 2, Wnorm_rev });
				}
				break;
			}
		}

		if (!_field_reverb_slice.Q_time.empty()) _field_reverb.Q_angle_time.insert({ ray.first, _field_reverb_slice });
	}

	return ret;
};

std::map<float, float> reverb::sum_spectrum(field_reverb&_field_reverb, float const &_t)
{
	std::map<float, float> specrtum;

	for (auto ray : _field_reverb.Q_angle_time)
	{
		for (auto step_ray : ray.second.Q_time)
		{
			int time = int(step_ray.first / _t) + 1;
			specrtum[_t * time] += step_ray.second;
		}
	}

// 	FILE *fff;
// 	fff = fopen("reverb.txt", "w");
// 	for (auto iter = specrtum.begin(); iter != specrtum.end(); iter++)
// 	{
// 		fprintf(fff, "%4.9f  %4.9f \n", iter->first, iter->second);
// 	}
// 	fclose(fff);

	return specrtum;
};
