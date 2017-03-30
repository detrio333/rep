#include "../headers/brdy.h"

HAC_USING_LIBRARY_NAMESPACE;

error const bdryMod::form_curve_bty_brdy(std::map<float, GDB::det::boundary> const &bottom_slices)
{
	error ret = error::OK_terminate;

	for (auto iter_r = bottom_slices.begin(); iter_r != bottom_slices.end(); iter_r++)
	{
		bdryMod::brdy _brdy;
		for (auto iter = iter_r->second.disp_refl.begin(); iter != iter_r->second.disp_refl.end(); iter++)
		{
			bdryMod::brdy::RefCoMod::ReflectionCoef _ReflectionCoef;
			_ReflectionCoef.R = pow(10.f, (iter->k_refl / 20.f));	// перевод из дЅ в отн.ед.
			_ReflectionCoef.phi = iter->phase_refl;
			_ReflectionCoef.theta = iter->angles;
			_brdy.RefCo.RefSet.insert({ _ReflectionCoef.theta, _ReflectionCoef });
			_brdy.ScaSet.insert({ iter->angles, pow(10.f, (iter->k_disp / 20.f)) });
		}
		_brdy.x = { iter_r->first, iter_r->second.depth_alt };
		_brdy.x = { iter_r->first, iter_r->second.depth_alt };
		map_brdy.insert({ _brdy.x.r, _brdy });
	}

	// учет возможых выходов за пределы расчетов
	{
		brdy brdy1 = (--map_brdy.end())->second;
		brdy1.x.r = sqrt(FLT_MAX);
		map_brdy.insert({ brdy1.x.r, brdy1 });

		brdy1 = map_brdy.begin()->second;
		brdy1.x.r = -sqrt(FLT_MAX);
		map_brdy.insert({ brdy1.x.r, brdy1 });
	}

	for (auto iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
	{
		coord_2d x = (++iter_brdy)->second.x;
		coord_2d y = (--iter_brdy)->second.x;
		iter_brdy->second.t = x - y;
		iter_brdy->second.Len = iter_brdy->second.t.abs();
		iter_brdy->second.t = iter_brdy->second.t / iter_brdy->second.Len;

		iter_brdy->second.normal.r = -iter_brdy->second.t.z;
		iter_brdy->second.normal.z = iter_brdy->second.t.r;
	}

	if (type == types::curve)	// curvilinear option : compute tangent and normal at node by averaging normals on adjacent segments
	{
		map_brdy.begin()->second.Nodet = { 1.0, 0.0 };	// tangent left - end  node
		(--map_brdy.end())->second.Nodet = { 1.0, 0.0 };// tangent right - end node

		for (auto iter_brdy = ++map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			float sss;
			//			sss = (--iter_brdy)->second.Len / (iter_brdy->second.Len + (++iter_brdy)->second.Len);
			sss = 0.5;	// странна вещь, но именнь так и написано
			iter_brdy->second.Nodet = (++iter_brdy)->second.t * sss + (--iter_brdy)->second.t * (1.0f - sss);
		}
		for (auto iter_brdy = map_brdy.begin(); iter_brdy != map_brdy.end(); iter_brdy++)
		{
			iter_brdy->second.Noden.r = -iter_brdy->second.Nodet.z;
			iter_brdy->second.Noden.z = iter_brdy->second.Nodet.r;
		}
		for (auto iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			// compute curvature in each segment
			float phi1 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy++;
			float phi2 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy--;
			iter_brdy->second.Kappa = (phi2 - phi1) / iter_brdy->second.Len;	// this is curvature = dphi / ds
		}
	}
	else
	{
		for (auto iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			iter_brdy->second.Kappa = 0;
		}
	}

	return ret;
};

error const bdryMod::form_curve_bty_brdy(const float &_depth, const float &_rBox)
{
	error ret = error::OK_terminate;

	type = types::linear;
	brdy _brdy;
	_brdy.x = { -_rBox, _depth };
	_brdy.t = { 1.f, 0.f };
	_brdy.normal = { 0.f, 1.f };
	_brdy.Kappa = 0;
	_brdy.Len = _brdy.t.abs();
	brdy::RefCoMod::ReflectionCoef _ReflectionCoef;
	_ReflectionCoef.R = 1;
	_ReflectionCoef.phi = static_cast<float>(M_PI);
	_ReflectionCoef.theta = 0;
	_brdy.RefCo.RefSet.insert({ 0.f, _ReflectionCoef });
	_ReflectionCoef.theta = static_cast<float>(M_PI_2);
	_brdy.RefCo.RefSet.insert({ static_cast<float>(M_PI_2), _ReflectionCoef });
	_brdy.ScaSet.insert({ 0.f, 0.05f });
	_brdy.ScaSet.insert({ static_cast<float>(M_PI_2), 0.05f });
	map_brdy.insert({ _brdy.x.r, _brdy });

	_brdy.x = { 0, _depth };
	map_brdy.insert({ _brdy.x.r, _brdy });

	_brdy.x = { _rBox, _depth };
	map_brdy.insert({ _brdy.x.r, _brdy });

	return ret;
};

error const bdryMod::form_curve_ati_brdy(std::map<float, GDB::det::boundary> const &surf_slices)
{
	error ret = error::OK_terminate;

	for (auto iter_r = surf_slices.begin(); iter_r != surf_slices.end(); iter_r++)
	{
		bdryMod::brdy _brdy;
		for (auto iter = iter_r->second.disp_refl.begin(); iter != iter_r->second.disp_refl.end(); iter++)
		{
			bdryMod::brdy::RefCoMod::ReflectionCoef _ReflectionCoef;
			_ReflectionCoef.R = pow(10.f, (iter->k_refl / 20.f));	// перевод из дЅ в отн.ед.
			_ReflectionCoef.phi = iter->phase_refl;
			_ReflectionCoef.theta = iter->angles;
			_brdy.RefCo.RefSet.insert({ _ReflectionCoef.theta, _ReflectionCoef });
			_brdy.ScaSet.insert({ iter->angles, pow(10.f, (iter->k_disp / 20.f)) });
		}
		_brdy.x = { iter_r->first, iter_r->second.depth_alt };
		map_brdy.insert({ _brdy.x.r, _brdy });
	}

	// учет возможых выходов за пределы расчетов
	{
		brdy brdy1 = (--map_brdy.end())->second;
		brdy1.x.r = sqrt(FLT_MAX);
		map_brdy.insert({ brdy1.x.r, brdy1 });

		brdy1 = map_brdy.begin()->second;
		brdy1.x.r = -sqrt(FLT_MAX);
		map_brdy.insert({ brdy1.x.r, brdy1 });
	}

	for (auto iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
	{
		coord_2d x = (++iter_brdy)->second.x;
		coord_2d y = (--iter_brdy)->second.x;
		iter_brdy->second.t = x - y;
		iter_brdy->second.Len = iter_brdy->second.t.abs();
		iter_brdy->second.t = iter_brdy->second.t / iter_brdy->second.Len;

		iter_brdy->second.normal.r = iter_brdy->second.t.z;
		iter_brdy->second.normal.z = -iter_brdy->second.t.r;
	}

	if (type == types::curve)	// curvilinear option : compute tangent and normal at node by averaging normals on adjacent segments
	{
		map_brdy.begin()->second.Nodet = { 1.0, 0.0 };	// tangent left - end  node
		(--map_brdy.end())->second.Nodet = { 1.0, 0.0 };// tangent right - end node

		for (auto iter_brdy = ++map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			float sss;
			//			sss = (--iter_brdy)->second.Len / (iter_brdy->second.Len + (++iter_brdy)->second.Len);
			sss = 0.5;	// странна вещь, но именнь так и написано
			iter_brdy->second.Nodet = (++iter_brdy)->second.t * sss + (--iter_brdy)->second.t * (1.0f - sss);
		}
		for (auto iter_brdy = map_brdy.begin(); iter_brdy != map_brdy.end(); iter_brdy++)
		{
			iter_brdy->second.Noden.r = iter_brdy->second.Nodet.z;
			iter_brdy->second.Noden.z = -iter_brdy->second.Nodet.r;
		}
		for (auto iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			// compute curvature in each segment
			float phi1 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy++;
			float phi2 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy--;
			iter_brdy->second.Kappa = (phi2 - phi1) / iter_brdy->second.Len;	// this is curvature = dphi / ds
		}
	}
	else
	{
		for (std::map<float, brdy>::iterator iter_brdy = map_brdy.begin(); iter_brdy != --map_brdy.end(); iter_brdy++)
		{
			iter_brdy->second.Kappa = 0;
		}
	}

	return ret;
};

error const bdryMod::form_curve_ati_brdy(const float &_depth, const float &_rBox, const float &_freq, const int &_heavy_sea)
{
	// коэфициенты отражени€ по умолчанию. в зависимости от волнени€ поверхности
	error ret = error::OK_terminate;

	type = types::linear;
	brdy _brdy;
	_brdy.x = { -_rBox, _depth };
	_brdy.t = { 1.f, 0.f };
	_brdy.normal = { 0.f, -1.f };
	_brdy.Kappa = 0;
	_brdy.Len = _brdy.t.abs();
	brdy::RefCoMod::ReflectionCoef _ReflectionCoef;
	auto d_alpha = static_cast<float>(M_PI) / 180;
	for (unsigned int i_alpha = 0; i_alpha < 91; i_alpha++)
	{
		_ReflectionCoef.theta = d_alpha*i_alpha;
		auto x = sqrt(static_cast<float>(M_PI)* 2);
		auto x1 = sqrt(x*_heavy_sea*_freq / 1000.f);
		auto x2 = pow(x*_heavy_sea, 0.1f);
		auto v2 = abs(1 - 0.56f*pow(x1, 3)*x2*sin(_ReflectionCoef.theta));
		auto pov = v2*v2;
		_ReflectionCoef.R = __min(pov, 1);
		_ReflectionCoef.phi = static_cast<float>(M_PI);
		_brdy.RefCo.RefSet.insert({ _ReflectionCoef.theta, _ReflectionCoef });
	}

	_brdy.ScaSet.insert({ 0.f, 0.05f });
	_brdy.ScaSet.insert({ static_cast<float>(M_PI_2), 0.05f });

	map_brdy.insert({ _brdy.x.r, _brdy });

	_brdy.x = { 0, _depth };
	map_brdy.insert({ _brdy.x.r, _brdy });

	_brdy.x = { _rBox, _depth };
	map_brdy.insert({ _brdy.x.r, _brdy });

	return ret;
};

bdryMod::brdy::RefCoMod::ReflectionCoef bdryMod::brdy::RefCoMod::REFCO(float &_alpha)
{
	ReflectionCoef RInt;

	auto iter_up = RefSet.upper_bound(_alpha);
	if (iter_up == RefSet.end())
	{
		iter_up--;
		RInt = iter_up->second;
		RInt.theta = _alpha;
	}
	auto iter_low = iter_up;
	iter_low--;

	auto Alpha = (_alpha - iter_low->second.theta) / (iter_up->second.theta - iter_low->second.theta);
	auto Alpha_m = 1 - Alpha;
	RInt.R = Alpha_m * iter_low->second.R + Alpha * iter_up->second.R;
	RInt.phi = Alpha_m * iter_low->second.phi + Alpha * iter_up->second.phi;
	RInt.theta = _alpha;

	return RInt;
};