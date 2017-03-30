#include "../headers/bellhop.h"

#include "math.h"

HAC_USING_LIBRARY_NAMESPACE;

error const bellhop::calc_ray_field()
{
	ss_at_source = local_hydrology->ssp.get_ssp(0.f, local_calc_param.depth_source);

	// забиваем начальную точку
	src_coord.r = 0;
	src_coord.z = local_calc_param.depth_source;

	// вычисление трасс лучей
	for (auto alpha_iter : src_beam.pow_in_angle)
	{
		auto ray = trace(alpha_iter);
		rays_field.insert({ alpha_iter, ray });
	}

	return error::OK_terminate;
};

error const bellhop::calc_clutter_field()
{
	std::map<float, std::vector<step_ray>> rays;
	ss_at_source = local_hydrology->ssp.get_ssp(0.f, local_calc_param.depth_source);

	// забиваем начальную точку
	src_coord.r = 0;
	src_coord.z = local_calc_param.depth_source;

	// вычисление трасс лучей
	for (auto alpha_iter : src_beam.pow_in_angle)
	{
		auto ray = clutter_trace(alpha_iter);
		rays.insert({ alpha_iter, ray });
	}

	return error::OK_terminate;
};

error const bellhop::get_level_field(unsigned int &_Nd, unsigned int &_Nr, std::ostream&_stream)
{
	error ret = error::OK_terminate;

	if (rays_field.empty())
	{
		ret = error::field_dont_calculated;
		return ret;
	}
	// сетка приемников
	SdRdR.Nrd = _Nd;
	SdRdR.rd = new float[SdRdR.Nrd];
	for (int i_rd = 0; i_rd < SdRdR.Nrd; i_rd++)
	{
		SdRdR.rd[i_rd] = (zBox - 100) / (SdRdR.Nrd - 1)*i_rd;
	}
	SdRdR.Nr = _Nr;
	SdRdR.r = new float[SdRdR.Nr];
	for (int i_r = 0; i_r < SdRdR.Nr; i_r++)
	{
		SdRdR.r[i_r] = rBox * 0.99f/ (SdRdR.Nr - 1)*i_r;
	}
	U = new std::complex<float>*[SdRdR.Nrd];
	for (int i_Nrd_per_range = 0; i_Nrd_per_range < SdRdR.Nrd; i_Nrd_per_range++)
	{
		U[i_Nrd_per_range] = new std::complex<float>[SdRdR.Nr];
	}

	switch (local_calc_param.type_beam)
	{
	case calc_param::types_beam::gauss:
		for (auto iter : rays_field)
		{
			InfluenceGeoGaussian(iter.second, iter.first);
		}
		break;
	case calc_param::types_beam::geom:
		for (auto iter : rays_field)
		{
			InfluenceGeoHat(iter.second, iter.first);
		}
		break;
	default:
		break;
	}

	scalep();
	for (auto i_z = 0; i_z < SdRdR.Nrd; i_z++)
	{
		for (auto i_r = 0; i_r < SdRdR.Nr; i_r++)
		{
			_stream << abs(U[i_z][i_r]) << " ";
		}
		_stream << "\n";
	}

	delete U;
	delete SdRdR.rd;
	delete SdRdR.r;

	return ret;
};

std::vector<step_ray> bellhop::trace(float _angle)
{
	std::vector<step_ray> ray;
	step_ray _step_ray;
	_step_ray.c = ss_at_source.c;	 // получить скорость звука
	_step_ray.x = src_coord;
	_step_ray.Tray.r = cos(_angle) / _step_ray.c;
	_step_ray.Tray.z = sin(_angle) / _step_ray.c;
	_step_ray.p.r = 1.0f;
	_step_ray.p.z = 0.0f;
	if (local_calc_param.type_beam == calc_param::types_beam::geom)
	{
		_step_ray.q.r = 0.0f;
		_step_ray.q.z = 0.0f;
	} 
	else
	{
		_step_ray.q.r = 0.0f;
		_step_ray.q.z = 1.0f;
	}
	_step_ray.tau = 0.0f;
	_step_ray.Rfa = 1;
	_step_ray.phase = 0.0f;
	_step_ray.Len = 0;
	_step_ray.event = step_ray::events::volume;

	{
		// identify the top segment above the source
		auto top_iter = local_hydrology->top_brdy.map_brdy.upper_bound(_step_ray.x.r);
		auto *top_u = &top_iter->second;
		top_iter--;
		auto *top_l = &top_iter->second;

		// identify the bottom segment below the source
		auto bot_iter = local_hydrology->bot_brdy.map_brdy.upper_bound(_step_ray.x.r);
		auto *bot_u = &bot_iter->second;
		bot_iter--;
		auto *bot_l = &bot_iter->second;

		// Trace the beam
		// (note that Reflect alters the step index)
		auto dEndTop = _step_ray.x - top_l->x;	// vector pointing from top    to ray
		auto dEndBot = _step_ray.x - bot_l->x;	// vector pointing from bottom to ray

		auto DistBegTop = dEndTop * top_l->normal;
		auto DistBegBot = dEndBot * bot_l->normal;

		//	!!!!note above distance is really a negative distance(throughout the code)

		if (DistBegTop >= 0 || DistBegBot >= 0)
		{
			return ray;	// source must be within the medium
		}

		for (auto i_step = 0; i_step < constants::MaxSteps; i_step++)
		{
			ray.push_back(_step_ray);

			_step_ray = step(_step_ray, *top_l, *top_u, *bot_l, *bot_u);
			
			if (_step_ray.x.r >= rBox || _step_ray.x.z >= zBox || _step_ray.x.r < 0)
			{
				break;
			}

			// New altimetry segment ?
			if (_step_ray.x.r < top_l->x.r || _step_ray.x.r > top_u->x.r)
			{
				top_iter = local_hydrology->top_brdy.map_brdy.upper_bound(_step_ray.x.r);
				top_u = &top_iter->second;
				top_iter--;
				top_l = &top_iter->second;
			}

			// New bathymetry  segment ?
			if (_step_ray.x.r < bot_l->x.r || _step_ray.x.r > bot_u->x.r)
			{
				bot_iter = local_hydrology->bot_brdy.map_brdy.upper_bound(_step_ray.x.r);
				bot_u = &bot_iter->second;
				bot_iter--;
				bot_l = &bot_iter->second;
			}

			//  *** Reflections ? ***
			// Tests that ray at step i is inside, and ray at step i + 1 is outside
			// to detect only a crossing from inside to outside

			dEndTop = _step_ray.x - top_l->x;	// vector pointing from top    to ray
			dEndBot = _step_ray.x - bot_l->x;	// vector pointing from bottom to ray

			auto DistEndTop = dEndTop * top_l->normal;
			auto DistEndBot = dEndBot * bot_l->normal;

			if ((DistBegTop > 0 && DistEndTop > 0) || (DistBegBot > 0 && DistEndBot > 0))
			{
				break;
			}

			if (DistBegTop < 0 && DistEndTop >= 0)	// test top reflection
			{
				coord_2d TopnInt, ToptInt;
				if (local_hydrology->top_brdy.type == bdryMod::types::curve)
				{
					auto sss = dEndTop * top_l->t / top_l->Len;	// proportional distance along segment
					TopnInt = top_l->Noden * (1 - sss) + top_u->Noden * sss;
					ToptInt = top_l->Nodet * (1 - sss) + top_u->Nodet * sss;
				}
				else
				{
					TopnInt = top_l->normal;	// normal is constant in a segment
					ToptInt = top_l->t;
				}

				ray.push_back(_step_ray);
				_step_ray = reflect_top(_step_ray, ToptInt, TopnInt, top_l->Kappa, top_l->RefCo);
				_step_ray.event = step_ray::events::surface;

				if (_step_ray.Rfa < 0.005)
				{
					break;
				}
			}
			else if (DistBegBot < 0 && DistEndBot >= 0)	// test bottom reflection
			{
				coord_2d BotnInt, BottInt;
				if (local_hydrology->bot_brdy.type == bdryMod::types::curve)
				{
					auto sss = dEndBot * bot_l->t / bot_l->Len;	// proportional distance along segment
					BotnInt = bot_l->Noden * (1 - sss) + bot_u->Noden * sss;
					BottInt = bot_l->Nodet * (1 - sss) + bot_u->Nodet * sss;
				}
				else
				{
					BotnInt = bot_l->normal;	// normal is constant in a segment
					BottInt = bot_l->t;
				}

				ray.push_back(_step_ray);
				_step_ray = reflect_bot(_step_ray, BottInt, BotnInt, bot_l->Kappa, bot_l->RefCo);
				_step_ray.event = step_ray::events::seabed;

				if (_step_ray.Rfa < 0.005)
				{
					break;
				}
			}

			DistBegTop = DistEndTop;
			DistBegBot = DistEndBot;
		}
	}

	return ray;
};

std::vector<step_ray> bellhop::clutter_trace(float _angle)
{
	std::vector<step_ray> ray;
	step_ray _step_ray;
	_step_ray.c = ss_at_source.c;	 // получить скорость звука
	_step_ray.x = src_coord;
	_step_ray.Tray.r = cos(_angle) / _step_ray.c;
	_step_ray.Tray.z = sin(_angle) / _step_ray.c;
	_step_ray.p.r = 1.0f;
	_step_ray.p.z = 0.0f;
	if (local_calc_param.type_beam == calc_param::types_beam::geom)
	{
		_step_ray.q.r = 0.0f;
		_step_ray.q.z = 0.0f;
	}
	else
	{
		_step_ray.q.r = 0.0f;
		_step_ray.q.z = 1.0f;
	}
	_step_ray.tau = 0.0f;
	_step_ray.Rfa = 1;
	_step_ray.phase = 0.0f;
	_step_ray.Len = 0;

	{
		// identify the top segment above the source
		auto top_iter = local_hydrology->top_brdy.map_brdy.upper_bound(_step_ray.x.r);
		auto top_u = &top_iter->second;
		top_iter--;
		auto top_l = &top_iter->second;

		// identify the bottom segment below the source
		auto bot_iter = local_hydrology->bot_brdy.map_brdy.upper_bound(_step_ray.x.r);
		auto bot_u = &bot_iter->second;
		bot_iter--;
		auto bot_l = &bot_iter->second;

		// Trace the beam
		// (note that Reflect alters the step index)
		auto dEndTop = _step_ray.x - top_l->x;	// vector pointing from top    to ray
		auto dEndBot = _step_ray.x - bot_l->x;	// vector pointing from bottom to ray

		auto DistBegTop = dEndTop * top_l->normal;
		auto DistBegBot = dEndBot * bot_l->normal;

		//	!!!!note above distance is really a negative distance(throughout the code)

		if (DistBegTop >= 0 || DistBegBot >= 0)
		{
			return ray;	// source must be within the medium
		}

		for (auto i_step = 0; i_step < constants::MaxSteps; i_step++)
		{
//			ray.push_back(_step_ray);

			_step_ray = step(_step_ray, *top_l, *top_u, *bot_l, *bot_u);

			if (_step_ray.x.r >= rBox || _step_ray.x.z >= zBox || _step_ray.x.r < 0)
			{
				break;
			}

			// New altimetry segment ?
			if (_step_ray.x.r < top_l->x.r || _step_ray.x.r > top_u->x.r)
			{
				top_iter = local_hydrology->top_brdy.map_brdy.upper_bound(_step_ray.x.r);
				top_u = &top_iter->second;
				top_iter--;
				top_l = &top_iter->second;
			}

			// New bathymetry  segment ?
			if (_step_ray.x.r < bot_l->x.r || _step_ray.x.r > bot_u->x.r)
			{
				bot_iter = local_hydrology->bot_brdy.map_brdy.upper_bound(_step_ray.x.r);
				bot_u = &bot_iter->second;
				bot_iter--;
				bot_l = &bot_iter->second;
			}

			//  *** Reflections ? ***
			// Tests that ray at step i is inside, and ray at step i + 1 is outside
			// to detect only a crossing from inside to outside

			dEndTop = _step_ray.x - top_l->x;	// vector pointing from top    to ray
			dEndBot = _step_ray.x - bot_l->x;	// vector pointing from bottom to ray

			auto DistEndTop = dEndTop * top_l->normal;
			auto DistEndBot = dEndBot * bot_l->normal;

			if ((DistBegTop > 0 && DistEndTop > 0) || (DistBegBot > 0 && DistEndBot > 0))
			{
				break;
			}

			if (DistBegTop < 0 && DistEndTop >= 0)	// test top reflection
			{
				coord_2d TopnInt, ToptInt;
				if (local_hydrology->top_brdy.type == bdryMod::types::curve)
				{
					auto sss = dEndTop * top_l->t / top_l->Len;	// proportional distance along segment
					TopnInt = top_l->Noden * (1 - sss) + top_u->Noden * sss;
					ToptInt = top_l->Nodet * (1 - sss) + top_u->Nodet * sss;
				}
				else
				{
					TopnInt = top_l->normal;	// normal is constant in a segment
					ToptInt = top_l->t;
				}

				_step_ray = reflect_top(_step_ray, ToptInt, TopnInt, top_l->Kappa, top_l->RefCo);
				ray.push_back(_step_ray);

				break;
			}
			else if (DistBegBot < 0 && DistEndBot >= 0)	// test bottom reflection
			{
				coord_2d BotnInt, BottInt;
				if (local_hydrology->bot_brdy.type == bdryMod::types::curve)
				{
					auto sss = dEndBot * bot_l->t / bot_l->Len;	// proportional distance along segment
					BotnInt = bot_l->Noden * (1 - sss) + bot_u->Noden * sss;
					BottInt = bot_l->Nodet * (1 - sss) + bot_u->Nodet * sss;
				}
				else
				{
					BotnInt = bot_l->normal;	// normal is constant in a segment
					BottInt = bot_l->t;
				}

//				ray.push_back(_step_ray);
				_step_ray = reflect_bot(_step_ray, BottInt, BotnInt, bot_l->Kappa, bot_l->RefCo);

				if (_step_ray.Rfa < 0.005)
				{
					break;
				}
			}

			DistBegTop = DistEndTop;
			DistBegBot = DistEndBot;
		}
	}

	return ray;
};

step_ray bellhop::step(const step_ray &_ray,
	const bdryMod::brdy &_xtop_lower_r, const bdryMod::brdy &_xtop_upper_r,
	const bdryMod::brdy &_xbot_lower_r, const bdryMod::brdy &_xbot_upper_r)
{
	step_ray ray1,ray2;

	// Does a single step along the ray
	// x denotes the ray coordinate, (r, z)
	// Tray denotes the scaled tangent to the ray

	// The numerical integrator used here is a version of the polygon(a.k.a.midpoint, leapfrog, or Box method), and similar
	// to the Heun(second order Runge - Kutta method).
	// However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

	// *** Phase 1 of modified polygon method(an Euler step) ***

	auto ssp_at_0 = local_hydrology->ssp.get_ssp(_ray.x.r, _ray.x.z);
	auto layer = local_hydrology->ssp.z_low;
	auto csq0 = ssp_at_0.c * ssp_at_0.c;
//	float cnn0_csq0 = ssp_at_0.crr * _ray.Tray.z * _ray.Tray.z - 2.f * ssp_at_0.crz * _ray.Tray.r * _ray.Tray.z + ssp_at_0.czz * _ray.Tray.r * _ray.Tray.r;

	auto h = deltas;   // initially set the step h, to the basic one, deltas
	reducestep(_ray.x, _ray.Tray, ssp_at_0.c, _xtop_lower_r, _xtop_upper_r, _xbot_lower_r, _xbot_upper_r, h);
	auto halfh = h * 0.5f;   // first step of the modified polygon method is a half step

	ray1.x = _ray.x + _ray.Tray * (halfh * ssp_at_0.c);
	ray1.Tray = _ray.Tray - ssp_at_0.grad_c * (halfh / csq0);
	ray1.p = _ray.p;// -_ray.q * (halfh * cnn0_csq0);
	ray1.q = _ray.q + _ray.p * (halfh * ssp_at_0.c);

	// *** Phase 2 of modified polygon method ***

	auto ssp_at_1 = local_hydrology->ssp.get_ssp(ray1.x.r, ray1.x.z);
	auto csq1 = ssp_at_1.c * ssp_at_1.c;
//	float cnn1_csq1 = ssp_at_1.crr * ray1.Tray.z * ray1.Tray.z - 2.f * ssp_at_1.crz * ray1.Tray.r * ray1.Tray.z + ssp_at_1.czz * ray1.Tray.r * ray1.Tray.r;

	if (h != min_step)
	{
		reducestep(_ray.x, ray1.Tray, ssp_at_1.c, _xtop_lower_r, _xtop_upper_r, _xbot_lower_r, _xbot_upper_r, h);
	}

	// use blend of f' based on proportion of a full step used.
	auto w1 = h / (2.0f * halfh);
	auto w0 = 1.0f - w1;
	auto hw0 = h * w0;
	auto hw1 = h * w1;

	auto a0 = hw0 * ssp_at_0.c;
	auto a1 = hw1 * ssp_at_1.c;

	ray2.x = _ray.x + _ray.Tray * (a0) + ray1.Tray * (a1);
	ray2.Tray = _ray.Tray - ssp_at_0.grad_c * (hw0 / csq0) - ssp_at_1.grad_c * (hw1 / csq1);
	ray2.p = _ray.p;// -_ray.q * (hw0 * cnn0_csq0) - ray1.q * (hw1 * cnn1_csq1);
	ray2.q = _ray.q + _ray.p * (a0) + ray1.p * (a1);
	ray2.tau = _ray.tau + hw0 / ssp_at_0.c + hw1 / ssp_at_1.c;
	ray2.Rfa = _ray.Rfa;
	ray2.phase = _ray.phase;
	auto tmp = _ray.x - ray2.x;
	ray2.Len = _ray.Len + tmp.abs();

	// If we crossed an interface, apply jump condition

	auto ssp_at_2 = local_hydrology->ssp.get_ssp(ray2.x.r, ray2.x.z);
	ray2.c = ssp_at_2.c;
	
	if (layer != local_hydrology->ssp.z_low)
	{
		coord_2d gradcjump = ssp_at_2.grad_c - ssp_at_0.grad_c;	// this is precise only for c - linear layers
		coord_2d ray2n(-ray2.Tray.z, ray2.Tray.r);

		auto cnjump = gradcjump * ray2n;
		auto csjump = gradcjump * ray2.Tray;

		auto RM = ray2.Tray.r / ray2.Tray.z;
		auto RN = -RM * (2 * cnjump - RM * csjump) / ray2.c;
		ray2.p = ray2.p + ray2.q * RN;
	}

	ray2.event = step_ray::events::volume;

	return ray2;
};

void bellhop::reducestep(const coord_2d &_ray_x, const coord_2d &_ray_tray, const float &_c,
	const bdryMod::brdy &_xtop_lower_r, const bdryMod::brdy &_xtop_upper_r,
	const bdryMod::brdy &_xbot_lower_r, const bdryMod::brdy &_xbot_upper_r,
	float & _h)	const //!< Reduces the ray step size to make sure we land on interfaces and boundaries
{
	// Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
	// Need to keep in mind possibility that user put source right on an interface
	// and that multiple events can occur(crossing interface, top, and bottom in a single step).
	// reminder: rayt is a scaled tangent; c * rayt is the unit tangent

	auto x = _ray_x + _ray_tray * (_h*_c);	// make a trial step

	// interface crossing in depth
	if (_ray_tray.z != 0)
	{
		float h1 = 1000;
		if (local_hydrology->ssp.z_low > x.z)
		{
			h1 = (local_hydrology->ssp.z_low - _ray_x.z) / (_c * _ray_tray.z);
		}
		else if (local_hydrology->ssp.z_up < x.z)
		{
			h1 = (local_hydrology->ssp.z_up - _ray_x.z) / (_c * _ray_tray.z);
		}
		if (h1 <= min_step)
		{
			_h = min_step;
			return;
		}

		_h = _h < h1 ? _h : h1;
	}


	// top crossing
	if ((x - _xtop_lower_r.x) * _xtop_lower_r.normal > 0.f)
	{
		auto e = _ray_x - _xtop_lower_r.x;
		auto h2 = -(e * _xtop_lower_r.normal) / (_c * ( _ray_tray * _xtop_lower_r.normal));
		_h = _h < h2 ? _h : h2;
	}

	// bottom crossing
	if ((x - _xbot_lower_r.x) * _xbot_lower_r.normal > 0.f)
	{
		auto e = _ray_x - _xbot_lower_r.x;
		auto h3 = -(e * _xbot_lower_r.normal) / (_c * (_ray_tray * _xbot_lower_r.normal));
		_h = _h < h3 ? _h : h3;
	}

	// top segment crossing in range
	if (_ray_tray.r != 0)
	{
		float h4;
		if (x.r < _xtop_lower_r.x.r)
		{
			h4 = -(_ray_x.r - _xtop_lower_r.x.r) / (_ray_tray.r*_c);
			_h = _h < h4 ? _h : h4;
		}
		else if (x.r > _xtop_upper_r.x.r)
		{
			h4 = -(_ray_x.r - _xtop_upper_r.x.r) / (_ray_tray.r*_c);
			_h = _h < h4 ? _h : h4;
		}
	}

	// bottom segment crossing in range
	if (_ray_tray.r != 0)
	{
		float h5;
		if (x.r < _xbot_lower_r.x.r)
		{
			h5 = -(_ray_x.r - _xbot_lower_r.x.r) / (_ray_tray.r*_c);
			_h = _h < h5 ? _h : h5;
		}
		else if (x.r > _xbot_upper_r.x.r)
		{
			h5 = -(_ray_x.r - _xbot_upper_r.x.r) / (_ray_tray.r*_c);
			_h = _h < h5 ? _h : h5;
		}
	}

	_h = _h > min_step ? _h : min_step;

};

step_ray bellhop::reflect_top(const step_ray  _ray, const coord_2d &_tbdry, const coord_2d &_nbdry, const float &_kappa, bdryMod::brdy::RefCoMod &_RefCo)
{
	step_ray ray_out;
	// here's the geometric part, changing the ray direction
	ray_out.x = _ray.x;
	ray_out.Len = _ray.Len;

	auto Tg = _tbdry * _ray.Tray;	// component of ray tangent, along boundary
	auto Th = _nbdry * _ray.Tray;	// component of ray tangent, normal to boundary

	ray_out.Tray = _ray.Tray - _nbdry * (2.f * Th);
	coord_2d rayn(-_ray.Tray.z, _ray.Tray.r);

	// Calculate the change in curvature
	// Based on formulas given by Muller, Geoph.J.R.A.S., 79 (1984).

	auto sound = local_hydrology->ssp.get_ssp(ray_out.x.r, ray_out.x.z);

	// rayn and ray%t should be multiplied by c to produce unit normals and tangents
	// Therefore cn, cs are off by that same factor

	auto cn = rayn * sound.grad_c;
	auto cs = sound.grad_c * _ray.Tray;

	auto RN = 2 * _kappa / sound.c / sound.c / Th;	// boundary curvature correction

	cn = -cn;	// flip sign for top reflection
	RN = -RN;

	auto RM = Tg / Th;
	RN = RN + RM * (4 * cn - 2 * RM * cs) / sound.c;	// dividing by c instead of c ^ 2 compensates for the missing factor in cn, cs

	ray_out.c = sound.c;
	ray_out.p = _ray.p + _ray.q * RN;
	ray_out.q = _ray.q;
	ray_out.tau = _ray.tau;
	// account for phase change

	// vacuum - для отладки
// 	ray_out.tau = _ray.tau;
// 	ray_out.Rfa = _ray.Rfa;
// 	ray_out.phase = _ray.phase + float(M_PI);

	// потом переделать на это 
	auto alpha = abs(atan2(Th, Tg));	// angle of incidence(relative to normal to bathymetry)
	if (alpha > static_cast<float>(M_PI_2))
	{
		alpha = static_cast<float>(M_PI)-alpha;	// reflection coefficient is symmetric about 90 degrees
	}
	auto RInt = _RefCo.REFCO(alpha);
	ray_out.Rfa = _ray.Rfa * RInt.R;
	ray_out.phase = _ray.phase + RInt.phi;

	return ray_out;
};

step_ray bellhop::reflect_bot(const step_ray _ray, const coord_2d &_tbdry, const coord_2d &_nbdry, const float &_kappa, bdryMod::brdy::RefCoMod &_RefCo)
{
	step_ray ray_out;
	// here's the geometric part, changing the ray direction
	ray_out.x = _ray.x;
	ray_out.Len = _ray.Len;

	auto Tg = _tbdry * _ray.Tray;	// component of ray tangent, along boundary
	auto Th = _nbdry * _ray.Tray;	// component of ray tangent, normal to boundary

	ray_out.Tray = _ray.Tray - _nbdry * (2.f * Th);
	coord_2d rayn(-_ray.Tray.z, _ray.Tray.r);

	// Calculate the change in curvature
	// Based on formulas given by Muller, Geoph.J.R.A.S., 79 (1984).

	auto sound = local_hydrology->ssp.get_ssp(ray_out.x.r, ray_out.x.z);

	// rayn and ray%t should be multiplied by c to produce unit normals and tangents
	// Therefore cn, cs are off by that same factor

	auto cn = rayn * sound.grad_c;
	auto cs = sound.grad_c * _ray.Tray;

	auto RN = 2 * _kappa / sound.c / sound.c / Th;	// boundary curvature correction

	auto RM = Tg / Th;
	RN = RN + RM * (4 * cn - 2 * RM * cs) / sound.c;	// dividing by c instead of c ^ 2 compensates for the missing factor in cn, cs

	ray_out.c = sound.c;
	ray_out.p = _ray.p + _ray.q * RN;
	ray_out.q = _ray.q;
	// account for phase change

	auto alpha = abs(atan2(Th, Tg));	// angle of incidence(relative to normal to bathymetry)
	if (alpha > static_cast<float>(M_PI_2))
	{
		alpha = static_cast<float>(M_PI)-alpha;	// reflection coefficient is symmetric about 90 degrees
	}
	auto RInt = _RefCo.REFCO(alpha);
	ray_out.Rfa = _ray.Rfa * RInt.R;
	ray_out.phase = _ray.phase + RInt.phi;

	return ray_out;
};

bellhop::bellhop()
{

};

bellhop::~bellhop()
{

};

error const bellhop::set_param(calc_param const &_local_calc_param)
{
	error ret = error::OK_terminate;

	if (local_hydrology == NULL)
	{
		ret = error::brdy_param_not_set;
		return ret;
	}

	local_calc_param = _local_calc_param;

 	rBox = local_calc_param.dist_max*1.01f;
	zBox = local_hydrology->bot_brdy.map_brdy.begin()->second.x.z;

	for (auto iter_r : local_hydrology->bot_brdy.map_brdy)
	{
		zBox = __max(zBox, iter_r.second.x.z);
	}
	zBox += 100;

	auto freq = local_hydrology->FB.get_freq_avg();
	switch (local_calc_param.type_attenuation)
	{
	case calc_param::types_attenuation::Thorpe:
	{
												  auto f2 = pow((freq / 1000), 2);
												  // Updated formula from JKPS Eq. 1.34
												  koef_attenuation = 0.0033f + 0.11f * f2 / (1.0f + f2) + 44.0f * f2 / (4100.0f + f2) + 0.00033f * f2;	// dB / km
												  koef_attenuation /= 8685.8896;	// Nepers / m
	}
		break;
	case calc_param::types_attenuation::database:
		koef_attenuation = local_hydrology->koef_attenuation / 8685.8896f;	// Nepers / m
		break;
	case calc_param::types_attenuation::none:
		koef_attenuation = 0;
	default:
		koef_attenuation = 0;
		break;
	}

	unsigned int Nbeams(__max(int(0.3f * freq * local_calc_param.dist_max / 1500), 50));	// automatically estimate NBeams to use
	if (!local_calc_param.src_beam.is_auto)
		Nbeams = local_calc_param.src_beam.pow_in_angle.size();	

	Dalpha = (local_calc_param.src_beam.angle.high - local_calc_param.src_beam.angle.low) / (Nbeams - 1);
	for (unsigned int i_alpha = 0; i_alpha < Nbeams; i_alpha++)
	{
		auto angle = local_calc_param.src_beam.angle.low + Dalpha*i_alpha;
		src_beam.pow_in_angle.push_back(angle);
	}
 
	deltas = zBox / 10;
 	min_step = deltas*float(1e-4) * 5;
 	
//	HSBot.BC = 'A';
//	HSBot.cP = { 1550.0f, 14.2006f };
//	HSBot.cS = { 0.0f, 0.0f };
//	HSBot.rho = 1.5f;

	return ret;
};

void bellhop::scalep()
{
	// Scale the pressure field

	//	if (RunType[0] /= 'C') U = SQRT(REAL(U)) !For incoherent run, convert intensity to pressure

	// add in attenuation
	std::complex<float> a(0, 0);
	for (auto ir = 0; ir < SdRdR.Nr; ir++)
	{
		float factor;
		if (local_calc_param.type_source == calc_param::types_source::cartesian)	// line source
		{
			factor = 4.0f * sqrt(static_cast<float>(M_PI));
		}
		else	// point source
		{
			if (SdRdR.r[ir] == 0)
			{
				factor = 0.f;	// avoid / 0 at origin, return pressure = 0
			}
			else
			{
				factor = -1.f*exp(-koef_attenuation*SdRdR.r[ir]) / sqrt(SdRdR.r[ir]);
			}
		}
		for (auto id = 0; id < SdRdR.Nrd; id++)
		{
			if (factor != 0 && U[id][ir] != a)
			{
				if (local_calc_param.type_coh != calc_param::types_coh::coherent)
				{
					U[id][ir] = sqrt(real(U[id][ir]))*factor;
				}
				else
				{
					U[id][ir] *= factor;
				}
			}
		}
	}
};

error const bellhop::get_ray_field(std::ostream&_stream) const
{
	error ret = error::OK_terminate;

	if (rays_field.empty())
	{
		ret = error::field_dont_calculated;
		return ret;
	}
	
	_stream.precision(5);
	for (auto iter1 : rays_field)
	{
		for (auto iter2 = iter1.second.begin(); iter2 != iter1.second.end(); iter2++)
		{
			_stream << iter2->x.r << " " << iter2->x.z << "\n";
		}
		for (auto iter2 = --iter1.second.end(); iter2 != iter1.second.begin(); iter2--)
		{
			_stream << iter2->x.r << " " << iter2->x.z << "\n";
		}
	}

	return ret;
};

error const bellhop::get_ray_field(std::list<std::map<float, float>> &_container) const
{
	error ret = error::OK_terminate;

	if (rays_field.empty())
	{
		ret = error::field_dont_calculated;
		return ret;
	}

	for (auto iter1 : rays_field)
	{
		std::map<float, float> map_coord;
		for (auto iter2 : iter1.second)
		{
			map_coord.insert({ iter2.x.r, iter2.x.z });
		}
		_container.push_back(map_coord);
	}

	return ret;
};

error const bellhop::get_gain(
	float& _r
	, float& _z
	, std::ostream& _stream
	) const		//!< функция расчёта передаточной характеристики поля для произвольного расположения приёмника
{
	error ret = error::OK_terminate;

	if (_r > rBox || _z > zBox || _r < 0 || _z < 0)
	{
		ret = error::limit_error;
		return ret;
	}

	if (rays_field.empty())
	{
		ret = error::field_dont_calculated;
		return ret;
	}

	switch (local_calc_param.type_beam)
	{
	case calc_param::types_beam::gauss:
		for (auto iter : rays_field)
		{
			InfluenceGeoGaussianOne(_stream, iter.second, iter.first, _r, _z);
		}
		break;
	case calc_param::types_beam::geom:
		for (auto iter : rays_field)
		{
			InfluenceGeoHatOne(_stream, iter.second, iter.first, _r, _z);
		}
		break;
	default:
		break;
	}

	return ret;
};

error const bellhop::get_gain(
	float& _r
	, float& _z
	, field_gain& _field_gain
	)		//!< функция расчёта передаточной характеристики поля для произвольного расположения приёмника
{
	error ret = error::OK_terminate;

	if (_r > rBox || _z > zBox || _r < 0 || _z < 0)
	{
		ret = error::limit_error;
		return ret;
	}

	if (rays_field.empty())
	{
		ret = error::field_dont_calculated;
		return ret;
	}

	_field_gain.c_em = ss_at_source.c;
	_field_gain.depth_em = local_calc_param.depth_source;
	_field_gain.depth_rcv = _z;
	_field_gain.dist = _r;
	_field_gain.freq = local_hydrology->FB.get_freq_avg();

	switch (local_calc_param.type_beam)
	{
	case calc_param::types_beam::gauss:
		for (auto iter : rays_field)
		{
			InfluenceGeoGaussianOne(_field_gain.rays, iter.second, iter.first, _r, _z);
		}
		break;
	case calc_param::types_beam::geom:
		for (auto iter : rays_field)
		{
			InfluenceGeoHatOne(_field_gain.rays, iter.second, iter.first, _r, _z);
		}
		break;
	default:
		break;
	}

	ss_at_point ss_at_receiver = local_hydrology->ssp.get_ssp(_r, _z);
	_field_gain.c_rcv = ss_at_receiver.c;

	return ret;
};

void bellhop::InfluenceGeoGaussianOne(std::multimap<float, field_gain::ray_desc> &_rays, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, Gaussian beams

	const int BeamWindow = 4;	// beam window : kills beams outside e^(-0.5 * ibwin^2)

	auto DS = sqrt(2.0f) * sin(omega * src_coord.z * _ray[0].Tray.z);	// Lloyd mirror pattern
	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	float Ratio1;
	// factor 1.2535 represents a sum of Gaussians in free space
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha))) / 1.2535f;	// point source
	}
	else
	{
		Ratio1 = 1 / 1.2535f;	// line  source
	}

	// add in attenuation
	float factor;
	if (local_calc_param.type_source == calc_param::types_source::cartesian)	// line source
	{
		factor = 4.0f * sqrt(static_cast<float>(M_PI));
	}
	else	// point source
	{
		if (_r == 0)
		{
			factor = 0.f;	// avoid / 0 at origin, return pressure = 0
		}
		else
		{
			factor = exp(-koef_attenuation*_r) / sqrt(_r);
		}
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto ray_1 = _ray[is];
		auto ray_0 = _ray[is - 1];

		auto rB = ray_1.x.r;

		// is r(ir) contained in[rA, rB}? Then compute beam influence
		if (_r < __min(rA, rB) || _r >= __max(rA, rB))
		{
			rA = rB;
			continue;
		}
		auto xray = ray_0.x;

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = ray_1.x - ray_0.x;
		auto rlen = rayt.abs();
		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = ray_0.q.r;
		if ((q <= 0.0f && qold > 0.0f) || (q >= 0.0f && qold < 0.0f)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto lambda = ray_0.c / (omega / (2 * static_cast<float>(M_PI)));
		auto sigma = __max(abs(ray_0.q.r), abs(ray_1.q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line
		auto sigma_l = __min(0.2f * (is + 1) * deltas / lambda, static_cast<float>(M_PI)* lambda);
		sigma = __max(sigma, sigma_l);
		auto RadMax = BeamWindow * sigma;

		if (RadMax > 200) continue;	// Савватеев, совсем мало мощности доходит при большом расплывании

		auto zmin = __min(ray_0.x.z, ray_1.x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(ray_0.x.z, ray_1.x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? 
		// If so, don't try to get depth limits: it's too complicated
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		coord_2d xrcvr(_r, _z);

		if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

		auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
		auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
		q = ray_0.q.r + s * (ray_1.q.r - ray_0.q.r);						// interpolated amplitude
		sigma = abs(q / q0);								// beam radius

		// посчитано это уже
		// calculate the beamwidth(must be at least pi * lambda, except in the nearfield)
		sigma = sigma > sigma_l ? sigma : sigma_l;//fmax(sigma, sigma_l);

		if (n < BeamWindow * sigma)	// Within beam window?
		{
			field_gain::ray_desc out_ray;
			auto A = abs(q0 / q);
			auto delay = ray_0.tau + s * (ray_1.tau - ray_0.tau);	// interpolated delay
			auto constt = Ratio1 * sqrt(ray_1.c / abs(q)) * ray_1.Rfa;
			auto phaseInt = phase;
			if ((q <= 0.0f && qold > 0.0f) || (q >= 0.0f && qold < 0.0f)) phaseInt = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
			if (local_calc_param.type_coh == calc_param::types_coh::semicoherent) constt *= DS;	// semi - coherent TL
			auto Amp = constt * exp(-0.5f * pow((n / sigma), 2)) / (2.f * sigma * A);

			out_ray.time = delay;
			out_ray.angle_em = _alpha;
			out_ray.angle_rcv = -atan2(ray_1.Tray.z, ray_1.Tray.r);
			out_ray.phase = ray_1.phase + phaseInt;
			out_ray.length = ray_1.Len;

			if (local_calc_param.type_coh == calc_param::types_coh::coherent)
			{
				std::complex<float> contri(0, -(omega * delay - ray_1.phase - phaseInt));	// фаза приходящего сигнала
				contri = exp(contri);
				out_ray.gain_coh = (Amp*factor)*contri;
				_rays.insert({ abs(out_ray.gain_coh), out_ray });
			}
			else // incoherent / semi - coherent TL
			{
				auto W = exp(-0.5f * pow((n / sigma), 2)) / (2.f * sigma * A);   //Gaussian decay
				out_ray.gain_incoh = Amp*factor / sqrt(W);
				_rays.insert({ out_ray.gain_incoh, out_ray });
			}
			break;
		}
		rA = rB;
	}
};

void bellhop::InfluenceGeoHatOne(std::multimap<float, field_gain::ray_desc> &_rays, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, hat - shaped beams

	auto DS = sqrt(2.0f) * sin(omega * src_coord.z * _ray[0].Tray.z);	// Lloyd mirror pattern
	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	float Ratio1;
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha)));	// point source
	}
	else
	{
		Ratio1 = 1;	// line  source
	}

	// add in attenuation
	float factor;
	if (local_calc_param.type_source == calc_param::types_source::cartesian)	// line source
	{
		factor = 4.0f * sqrt(static_cast<float>(M_PI));
	}
	else	// point source
	{
		if (_r == 0)
		{
			factor = 0.f;	// avoid / 0 at origin, return pressure = 0
		}
		else
		{
			factor = exp(-koef_attenuation*_r) / sqrt(_r);
		}
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto ray_1 = _ray[is];
		auto ray_0 = _ray[is - 1];

		auto rB = ray_1.x.r;

		if (_r < __min(rA, rB) || _r >= __max(rA, rB))
		{
			rA = rB;
			continue;
		}

		auto xray = ray_0.x;
		if (abs(rB - rA) < FLT_EPSILON) continue;	// jump to next step if duplicate point

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = ray_1.x - ray_0.x;
		auto rlen = rayt.abs();
		// if (rlen < deltas*1e-4 * 5) continue;	// if duplicate point in ray, skip to next step along the ray
		if (rlen < 100 * FLT_EPSILON) continue;	// if duplicate point in ray, skip to next step along the ray

		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = ray_0.q.r;
		if ((q <= 0.0 && qold > 0.0) || (q >= 0.0 && qold < 0.0)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto RadMax = __max(abs(ray_0.q.r), abs(ray_1.q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line

		if (RadMax > 100) continue;	// Савватеев, совсем мало мощности доходит при большом расплывании

		auto zmin = __min(ray_0.x.z, ray_1.x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(ray_0.x.z, ray_1.x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? Then don't try to get depth limits: it's too complicated
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		// computed beam influence for this segment of the ray
		coord_2d xrcvr(_r, _z);

		if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

		auto dqds = ray_1.q.r - ray_0.q.r;
		auto dtauds = ray_1.tau - ray_0.tau;

		auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
		auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
		q = ray_0.q.r + s * dqds;						// interpolated amplitude
		RadMax = abs(q / q0);								// beam radius
		//		RadMax = __max(RadMax, 100);

		if (n < RadMax)
		{
			field_gain::ray_desc out_ray;
			auto A = 1 / RadMax;
			auto delay = ray_0.tau + s * dtauds;	// interpolated delay
			auto constt = Ratio1 * sqrt(ray_1.c / abs(q)) * A * ray_1.Rfa;
			auto phaseInt = phase;
			if ((q <= 0.0 && qold > 0.0) || (q >= 0.0 && qold < 0.0)) phaseInt = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics

			if (local_calc_param.type_coh == calc_param::types_coh::semicoherent) constt *= DS;	// semi - coherent TL
			auto Amp = constt * (RadMax - n);

			out_ray.time = delay;
			out_ray.angle_em = _alpha;
			out_ray.angle_rcv = -atan2(ray_1.Tray.z, ray_1.Tray.r);
			out_ray.phase = ray_1.phase + phaseInt;
			out_ray.length = ray_1.Len;

			if (local_calc_param.type_coh == calc_param::types_coh::coherent)
			{
				std::complex<float> contri(0, -(omega * delay - _ray[is].phase - phaseInt));	// фаза приходящего сигнала
				contri = exp(contri);
				out_ray.gain_coh = (Amp*factor)*contri;
				_rays.insert({ abs(out_ray.gain_coh), out_ray });
			}
			else // incoherent / semi - coherent TL
			{
				auto W = (RadMax - n) / RadMax;	// hat function : 1 on center, 0 on edge
				out_ray.gain_incoh = Amp*factor / sqrt(W);
				_rays.insert({ out_ray.gain_incoh, out_ray });
			}
			break;
		}
		rA = rB;
	}
};

void bellhop::InfluenceGeoGaussianOne(std::ostream& _stream, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, Gaussian beams

	const int BeamWindow = 4;	// beam window : kills beams outside e^(-0.5 * ibwin^2)

	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	float Ratio1;
	// factor 1.2535 represents a sum of Gaussians in free space
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha))) / 1.2535f;	// point source
	}
	else
	{
		Ratio1 = 1 / 1.2535f;	// line  source
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto ray_1 = _ray[is];
		auto ray_0 = _ray[is - 1];

		float rB = ray_1.x.r;

		// is r(ir) contained in[rA, rB}? Then compute beam influence
		if (_r < __min(rA, rB) || _r >= __max(rA, rB))
		{
			rA = rB;
			continue;
		}
		auto xray = ray_0.x;

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = ray_1.x - ray_0.x;
		auto rlen = rayt.abs();
		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = ray_0.q.r;
		if ((q <= 0.0f && qold > 0.0f) || (q >= 0.0f && qold < 0.0f)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto lambda = ray_0.c / (omega / (2 * static_cast<float>(M_PI)));
		auto sigma = __max(abs(ray_0.q.r), abs(ray_1.q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line
		auto sigma_l = __min(0.2f * (is + 1) * deltas / lambda, static_cast<float>(M_PI)* lambda);
		sigma = __max(sigma, sigma_l);
		auto RadMax = BeamWindow * sigma;

		if (RadMax > 200) continue;	// Савватеев, совсем мало мощности доходит при большом расплывании

		auto zmin = __min(ray_0.x.z, ray_1.x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(ray_0.x.z, ray_1.x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? 
		// If so, don't try to get depth limits: it's too complicated
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		coord_2d xrcvr(_r, _z);

		if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

		auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
		auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
		q = ray_0.q.r + s * (ray_1.q.r - ray_0.q.r);						// interpolated amplitude
		sigma = abs(q / q0);								// beam radius

		// посчитано это уже
		// calculate the beamwidth(must be at least pi * lambda, except in the nearfield)
		sigma = sigma > sigma_l ? sigma : sigma_l;//fmax(sigma, sigma_l);

		if (n < BeamWindow * sigma)	// Within beam window?
		{
			for (auto iter2 = _ray.begin(); iter2 != _ray.end(); iter2++)
			{
				_stream << iter2->x.r << " " << iter2->x.z << "\n";
			}
			for (auto iter2 = --_ray.end(); iter2 != _ray.begin(); iter2--)
			{
				_stream << iter2->x.r << " " << iter2->x.z << "\n";
			}

			break;
		}
		rA = rB;
	}
};

void bellhop::InfluenceGeoHatOne(std::ostream& _stream, const std::vector<step_ray> &_ray, const float &_alpha, const float&_r, const float&_z) const
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, hat - shaped beams

	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	float Ratio1;
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha)));	// point source
	}
	else
	{
		Ratio1 = 1;	// line  source
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto ray_1 = _ray[is];
		auto ray_0 = _ray[is - 1];

		auto rB = ray_1.x.r;

		if (_r < __min(rA, rB) || _r >= __max(rA, rB))
		{
			rA = rB;
			continue;
		}

		auto xray = ray_0.x;
		if (abs(rB - rA) < FLT_EPSILON) continue;	// jump to next step if duplicate point

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = ray_1.x - ray_0.x;
		auto rlen = rayt.abs();
		// if (rlen < deltas*1e-4 * 5) continue;	// if duplicate point in ray, skip to next step along the ray
		if (rlen < 100 * FLT_EPSILON) continue;	// if duplicate point in ray, skip to next step along the ray

		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = ray_0.q.r;
		if ((q <= 0.0 && qold > 0.0) || (q >= 0.0 && qold < 0.0)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto RadMax = __max(abs(ray_0.q.r), abs(ray_1.q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line

		if (RadMax > 100) continue;	// Савватеев, совсем мало мощности доходит при большом расплывании

		auto zmin = __min(ray_0.x.z, ray_1.x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(ray_0.x.z, ray_1.x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? Then don't try to get depth limits: it's too complicated
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		// computed beam influence for this segment of the ray
		coord_2d xrcvr(_r, _z);

		if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

		auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
		auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
		q = ray_0.q.r + s * (ray_1.q.r - ray_0.q.r);						// interpolated amplitude
		RadMax = abs(q / q0);								// beam radius

//		_stream.precision(4);
		if (n < RadMax)
		{
			for (auto iter2 = _ray.begin(); iter2 != _ray.end(); iter2++)
			{
				_stream << iter2->x.r << " " << iter2->x.z << "\n";
			}
			for (auto iter2 = --_ray.end(); iter2 != _ray.begin(); iter2--)
			{
				_stream << iter2->x.r << " " << iter2->x.z << "\n";
			}

			break;
		}
		rA = rB;
	}
};

void bellhop::InfluenceGeoGaussian(const std::vector<step_ray> &_ray, const float &_alpha)
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, Gaussian beams

	const int BeamWindow = 4;	// beam window : kills beams outside e^(-0.5 * ibwin^2)

	auto DS = sqrt(2.0f) * sin(omega * src_coord.z * _ray[0].Tray.z);	// Lloyd mirror pattern
	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	// what if never satistified ?
	int ir = 0;
	int irTT;
	for (int i_r = 0; i_r < SdRdR.Nr; i_r++)	// find index of first receiver to the right of rA
	{
		if (SdRdR.r[i_r] > rA)
		{
			ir = i_r;
			break;
		}
	}
	if (_ray[0].Tray.r < 0.f)	// if ray is traveling to the left, then we want the first receiver to the left of rA
	{
		ir--;
	}

	float Ratio1;
	// factor 1.2535 represents a sum of Gaussians in free space
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha))) / 1.2535f;	// point source
	}
	else
	{
		Ratio1 = 1 / 1.2535f;	// line  source
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto ray_1 = _ray[is];
		auto ray_0 = _ray[is - 1];
		auto rB = ray_1.x.r;
		auto xray = ray_0.x;

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = ray_1.x - ray_0.x;
		auto rlen = rayt.abs();
		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = ray_0.q.r;
		if ((q <= 0.0f && qold > 0.0f) || (q >= 0.0f && qold < 0.0f)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto lambda = ray_0.c / (omega / (2 * static_cast<float>(M_PI)));
		auto sigma = __max(abs(ray_0.q.r), abs(ray_1.q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line
		auto sigma_l = __min(0.2f * (is + 1) * deltas / lambda, static_cast<float>(M_PI)* lambda);
		sigma = __max(sigma, sigma_l);
		auto RadMax = BeamWindow * sigma;

		auto zmin = __min(ray_0.x.z, ray_1.x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(ray_0.x.z, ray_1.x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? 
		// If so, don't try to get depth limits: it's too complicated
		// крутой луч?
		// Если да, то не пытайся получить ограничение по глубине: слишком сложно
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		auto dqds = ray_1.q.r - ray_0.q.r;
		auto dtauds = ray_1.tau - ray_0.tau;

		// computed beam influence for this segment of the ray
		while (true)
		{
			// is r(ir) contained in[rA, rB}? Then compute beam influence
			if (SdRdR.r[ir] >= __min(rA, rB) && SdRdR.r[ir] < __max(rA, rB))
			{
				for (int id = 0; id < SdRdR.Nrd; id++)	// Loop over receiver depths
				{
					coord_2d xrcvr = { SdRdR.r[ir], SdRdR.rd[id] };	// irregular   grid

					if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

					auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
					auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
					q = ray_0.q.r + s * dqds;						// interpolated amplitude
					sigma = abs(q / q0);								// beam radius

					// посчитано это уже
					// 					// calculate the beamwidth(must be at least pi * lambda, except in the nearfield)
					// 					lambda = _ray[is - 1].c / (omega / (2 * float(M_PI)));
					sigma = sigma > sigma_l ? sigma : sigma_l;//fmax(sigma, sigma_l);

					if (n < BeamWindow * sigma)	// Within beam window?
					{
						auto A = abs(q0 / q);
						auto delay = ray_0.tau + s * dtauds;	// interpolated delay
						auto constt = Ratio1 * sqrt(ray_1.c / abs(q)) * ray_1.Rfa;
						auto phaseInt = phase;
						if ((q <= 0.0f && qold > 0.0f) || (q >= 0.0f && qold < 0.0f)) phaseInt = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
						if (local_calc_param.type_coh == calc_param::types_coh::semicoherent) constt *= DS;	// semi - coherent TL
						auto Amp = constt * exp(-0.5f * pow((n / sigma), 2)) / (2.f * sigma * A);

						if (local_calc_param.type_coh == calc_param::types_coh::coherent)
						{
							std::complex<float> contri(0, -(omega * delay - ray_1.phase - phaseInt));
							contri = exp(contri);
							U[id][ir] += Amp*contri;

						}
						else // incoherent / semi - coherent TL
						{
							auto W = exp(-0.5f * pow((n / sigma), 2)) / (2.f * sigma * A);   //Gaussian decay
							U[id][ir] += pow((Amp / W), 2)*W;
						}
					}
				}
				// Next receiver depth
			}

			// bump receiver index, ir, towards rB
			if (SdRdR.r[ir] < rB)
			{
				if (ir >= SdRdR.Nr - 1) break;	// jump out of the search and go to next step on ray
				irTT = ir + 1;	// bump right
				if (SdRdR.r[irTT] >= rB) break;
			}
			else
			{
				if (ir <= 0) break;	// jump out of the search and go to next step on ray
				irTT = ir - 1;	// bump left
				if (SdRdR.r[irTT] <= rB) break;
			}
			ir = irTT;
			// Next receiver range
		}
		rA = rB;
	}
};

void bellhop::InfluenceGeoHat(const std::vector<step_ray> &_ray, const float &_alpha)
{
	// Computes the beam influence, i.e.
	// the contribution of a single beam to the complex pressure
	// This version uses geometric, hat - shaped beams

	auto DS = sqrt(2.0f) * sin(omega * src_coord.z * _ray[0].Tray.z);	// Lloyd mirror pattern
	auto q0 = _ray[0].c / Dalpha;	// Reference for J = q0 / q
	auto phase = 0.0f;
	auto qold = _ray[0].q.r;
	auto rA = _ray[0].x.r;	// range at start of ray

	int ir = 0;
	int irTT;
	for (int i_r = 0; i_r < SdRdR.Nr; i_r++)	// find index of first receiver to the right of rA
	{
		if (SdRdR.r[i_r] > rA)
		{
			ir = i_r;
			break;
		}
	}
	if (_ray[0].Tray.r < 0.f)	// if ray is traveling to the left, then we want the first receiver to the left of rA
	{
		ir--;
	}

	float Ratio1;
	if (local_calc_param.type_source == calc_param::types_source::cylindrical)
	{
		Ratio1 = sqrt(abs(cos(_alpha)));	// point source
	}
	else
	{
		Ratio1 = 1;	// line  source
	}

	for (unsigned int is = 1; is < _ray.size(); is++)	// Loop over steps
	{
		auto rB = _ray[is].x.r;
		auto xray = _ray[is - 1].x;
		if (abs(rB - rA) < FLT_EPSILON) continue;	// jump to next step if duplicate point

		// initialize the index of the receiver range
		if (is == 1)
		{
			if (rB > rA)	// ray is moving right
			{
				ir = 0;	// index all the way to the left
			}
			else	// ray is moving left
			{
				ir = SdRdR.Nr - 1;	// index all the way to the right
			}
		}

		// compute normalized tangent(compute it because we need to measure the step length)
		auto rayt = _ray[is].x - _ray[is - 1].x;
		auto rlen = rayt.abs();
		// if (rlen < deltas*1e-4 * 5) continue;	// if duplicate point in ray, skip to next step along the ray
		if (rlen < 100 * FLT_EPSILON) continue;	// if duplicate point in ray, skip to next step along the ray

		rayt = rayt / rlen;	// unit tangent to ray
		coord_2d rayn(-rayt.z, rayt.r);	// unit normal  to ray

		// phase shifts at caustics
		auto q = _ray[is - 1].q.r;
		if ((q <= 0.0 && qold > 0.0) || (q >= 0.0 && qold < 0.0)) phase = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics
		qold = q;

		auto RadMax = __max(abs(_ray[is - 1].q.r), abs(_ray[is].q.r)) / q0 / abs(rayt.r);	// beam radius projected onto vertical line
		auto zmin = __min(_ray[is - 1].x.z, _ray[is].x.z) - RadMax;	// min depth of ray segment
		auto zmax = __max(_ray[is - 1].x.z, _ray[is].x.z) + RadMax;	// max depth of ray segment

		// is this a steep ray ? Then don't try to get depth limits: it's too complicated
		if (abs(rayt.r) < 0.5)
		{
			zmin = -FLT_MAX;
			zmax = +FLT_MAX;
		}

		auto dqds = _ray[is].q.r - _ray[is - 1].q.r;
		auto dtauds = _ray[is].tau - _ray[is - 1].tau;

		// computed beam influence for this segment of the ray
		while (true)
		{
			// is r(ir) contained in[rA, rB}? Then compute beam influence
			// IF(ABS(r(ir) - rA) + ABS(rB - r(ir)) <= ABS(rB - rA)) THEN
			if (SdRdR.r[ir] >= __min(rA, rB) && SdRdR.r[ir] < __max(rA, rB))
			{
				for (int id = 0; id < SdRdR.Nrd; id++)	// Loop over receiver depths
				{
					coord_2d xrcvr(SdRdR.r[ir], SdRdR.rd[id]);

					if (xrcvr.z < zmin || xrcvr.z > zmax) continue;

					auto s = rayt * (xrcvr - xray) / rlen;	// proportional distance along ray
					auto n = abs(rayn * (xrcvr - xray));	// normal distance to ray
					q = _ray[is - 1].q.r + s * dqds;						// interpolated amplitude
					RadMax = abs(q / q0);								// beam radius

					if (n < RadMax)
					{
						auto A = 1 / RadMax;
						auto delay = _ray[is - 1].tau + s * dtauds;		// interpolated delay
						auto constt = Ratio1 * sqrt(_ray[is].c / abs(q)) * A * _ray[is].Rfa;
						auto phaseInt = phase;
						if ((q <= 0.0 && qold > 0.0) || (q >= 0.0 && qold < 0.0)) phaseInt = phase + static_cast<float>(M_PI_2);	// phase shifts at caustics

						if (local_calc_param.type_coh == calc_param::types_coh::semicoherent) constt *= DS;	// semi - coherent TL
						auto Amp = constt * (RadMax - n);

						if (local_calc_param.type_coh == calc_param::types_coh::coherent)
						{
							std::complex<float> contri(0, -(omega * delay - _ray[is].phase - phaseInt));
							contri = exp(contri);
							U[id][ir] += Amp*contri;

						}
						else // incoherent / semi - coherent TL
						{
							auto W = (RadMax - n) / RadMax;	// hat function : 1 on center, 0 on edge
							U[id][ir] += pow((Amp / W), 2)*W;
						}
					}
				}
				// Next receiver depth
			}

			// bump receiver index, ir, towards rB
			if (SdRdR.r[ir] < rB)
			{
				if (ir >= SdRdR.Nr - 1)	// jump out of the search and go to next step on ray
				{
					break;
				}
				irTT = ir + 1;	// bump right
				if (SdRdR.r[irTT] >= rB)
				{
					break;
				}
			}
			else
			{
				if (ir <= 0)	// jump out of the search and go to next step on ray
				{
					break;
				}
				irTT = ir - 1;	// bump left
				if (SdRdR.r[irTT] <= rB)
				{
					break;
				}
			}
			ir = irTT;
			// Next receiver range
		}
		rA = rB;
	}
};
