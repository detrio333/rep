#include "../headers/clutter.h"

HAC_USING_LIBRARY_NAMESPACE;

error const clutter::set_param(calc_param const &_local_calc_param)
{
	error ret = error::OK_terminate;

	if (local_hydrology == NULL)
	{
		ret = error::brdy_param_not_set;
		return ret;
	}

	local_calc_param = _local_calc_param;

	return ret;
};

error const clutter::get_clutter(
	const float &_depth_receiver,		//!< глубина приёмника, м
	field_clutter &local_clutter
	)const
{
	unsigned int Nj;	// количество испускаемых лучей

	switch (local_calc_param.clutter_type)
	{
	case calc_param::clutter_types::hard:
		Nj = 200;
		break;
	case calc_param::clutter_types::normal:
		Nj = 100;
		break;
	case calc_param::clutter_types::low:
		Nj = 50;
		break;
	default:
		Nj = 100;
		break;
	}

	float max_ugol = static_cast<float>(M_PI_2)-0.1f;

	bellhop _bellhop;
	// надо перевернуть потом
	_bellhop.local_hydrology = local_hydrology;

	bellhop::calc_param trace_calc_param;
	trace_calc_param.depth_source = _depth_receiver;
	trace_calc_param.dist_max = 50 * 1000;

	trace_calc_param.src_beam.pow_in_angle.resize(Nj);
	trace_calc_param.src_beam.is_auto = false;
	trace_calc_param.src_beam.angle.low = -max_ugol;
	trace_calc_param.src_beam.angle.high = max_ugol;


	_bellhop.set_param(trace_calc_param);

	//	Трассировка лучей
	_bellhop.calc_clutter_field();
	auto rays = _bellhop.rays_field;

	auto koef_attenuation = local_hydrology->koef_attenuation / 8685.8896f;
	decltype(koef_attenuation) com = 0;
	for (auto ray_iter = rays.begin(); ray_iter != rays.end(); ray_iter++)
	{
		if (ray_iter->second.empty())
		{
			local_clutter.Q.insert({ -ray_iter->first, 0.f });
		}
		else
		{
			step_ray _step_ray = *ray_iter->second.begin();
			float power = _step_ray.Rfa * _step_ray.Rfa * exp(-2 * _step_ray.Len*koef_attenuation) * _step_ray.Tray.z * _step_ray.c;
			//float power = _step_ray.Rfa * _step_ray.Rfa * pow(10.f, -0.1f * local_hydrology->koef_attenuation*_step_ray.x.r / 1000) * _step_ray.Tray.z * _step_ray.c;
			local_clutter.Q.insert({ -ray_iter->first, power });
			com += power;
		}
	}

	for (auto iter : local_clutter.Q)
	{
		iter.second /= com;
	}

	// взяты значения из луны
	switch (local_hydrology->heavy_sea)
	{
	case 0:	local_clutter.P = 9.25f;	break;
	case 1:	local_clutter.P = 17.5f;	break;
	case 2:	local_clutter.P = 29.54f;	break;
	case 3:	local_clutter.P = 33.06f;	break;
	case 4:	local_clutter.P = 35.56f;	break;
	case 5:	local_clutter.P = 38.6f;	break;
	case 6:	local_clutter.P = 41.2f;	break;
	case 7:	local_clutter.P = 44.08f;	break;
	case 8:	local_clutter.P = 45.8f;	break;
	case 9:	local_clutter.P = 50.6f;	break;
	default:break;
	}
	// добавим зависимость от глубины. хотя это всё условно
	if (_depth_receiver < 500)
	{
		local_clutter.P *= (1 - 0.1f*_depth_receiver / 500);
	}
	else
	{
		local_clutter.P *= 0.9f;
	}

	local_clutter.P = 4 * pow(10.f, local_clutter.P / 10 - 10);	// перевод дБ в Па^2

	return error::OK_terminate;
};