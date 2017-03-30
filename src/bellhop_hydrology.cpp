#include "../headers/bellhop_hydrology.h"

HAC_USING_LIBRARY_NAMESPACE;

error const bellhop_hydrology::set_param(const GDB::det::hydrology &_hydrology, const bdryMod::types &_ati_type, const bdryMod::types &_bty_type)
{
	error ret = error::OK_terminate;

	heavy_sea = _hydrology.W;

	for (auto volume : _hydrology.vol_slices)
	{
		decltype(scat.begin()->second) tmp_map;
		for (auto tmp : volume.second.scat)
		{
			tmp_map.insert({ tmp.first, pow(10.f, (tmp.second / 20.f)) });
		}
		auto iter = tmp_map.begin();
		tmp_map.insert({ iter->first - 10, iter->second });
		iter = --tmp_map.end();
		tmp_map.insert({ iter->first * 1.1f, iter->second });
		scat.insert({ volume.first, tmp_map });
	}

	FB = { _hydrology.freq, _hydrology.freq };

	if (_hydrology.bottom_slices.empty())
	{
		return error::brdy_param_not_set;
	}

	if (_hydrology.vol_slices.empty())
	{
		koef_attenuation = 0.f;
	} 
	else
	{
		koef_attenuation = _hydrology.vol_slices.begin()->second.att;	// пока используем только одно значение
	}

	// ВРСЗ
	ssp.set_ssp(_hydrology.ssp_slices);

	bot_brdy.type = _bty_type;
	bot_brdy.form_curve_bty_brdy(_hydrology.bottom_slices);

	if (_hydrology.surf_slices.empty())
	{
		top_brdy.form_curve_ati_brdy(0.f, (--_hydrology.bottom_slices.end())->first*1.01f, FB.get_freq_avg(), heavy_sea);
	}
	else
	{
		top_brdy.type = _ati_type;
		top_brdy.form_curve_ati_brdy(_hydrology.surf_slices);
	}

	return ret;
};

bellhop_hydrology bellhop_hydrology::reverse(const float &_r)
{
	bellhop_hydrology bellhop_hydrology_out;

	bellhop_hydrology_out.FB = this->FB;
	bellhop_hydrology_out.koef_attenuation = this->koef_attenuation;

	// врсз
	bellhop_hydrology_out.ssp.dist_number = this->ssp.dist_number;
	if (this->ssp.dist_number == sspmod::dist_numbers::one)
	{
		for (auto depth_ssp : this->ssp.depth_ssp)
		{
			bellhop_hydrology_out.ssp.depth_ssp.insert({ depth_ssp.first, depth_ssp.second });
		}
	}
	else
	{
		for (auto iter_r : this->ssp.range_depth_ssp)
		{
			if (iter_r.first < _r)
			{
				auto r = _r - iter_r.first;
				decltype(bellhop_hydrology_out.ssp.range_depth_ssp.begin()->second) depth_ssp;
				for (auto iter_z : iter_r.second)
				{
					decltype(depth_ssp.begin()->second) _sonic;
					_sonic = iter_z.second;
					_sonic.r = r;
					depth_ssp.insert({ _sonic.h, _sonic });
				}
				bellhop_hydrology_out.ssp.range_depth_ssp.insert({ r, depth_ssp });
			}
		}
	}

	// поверхность
	bellhop_hydrology_out.top_brdy.type = this->top_brdy.type;
	for (auto iter_r : this->top_brdy.map_brdy)
	{
		auto _brdy = iter_r.second;
		_brdy.x.r = _r - iter_r.first;

		bellhop_hydrology_out.top_brdy.map_brdy.insert({ _brdy.x.r, _brdy });
		if (iter_r.first > _r)
		{
			_brdy.x.r = 0;

			bellhop_hydrology_out.top_brdy.map_brdy.insert({ _brdy.x.r, _brdy });
			break;
		}
	}

	for (auto iter_brdy = this->top_brdy.map_brdy.begin(); iter_brdy != --this->top_brdy.map_brdy.end(); iter_brdy++)
	{
		coord_2d x = (++iter_brdy)->second.x;
		coord_2d y = (--iter_brdy)->second.x;
		iter_brdy->second.t = x - y;
		iter_brdy->second.Len = iter_brdy->second.t.abs();
		iter_brdy->second.t = iter_brdy->second.t / iter_brdy->second.Len;

		iter_brdy->second.normal.r = iter_brdy->second.t.z;
		iter_brdy->second.normal.z = -iter_brdy->second.t.r;
	}

	if (bellhop_hydrology_out.top_brdy.type == bdryMod::types::curve)	// curvilinear option : compute tangent and normal at node by averaging normals on adjacent segments
	{
		bellhop_hydrology_out.top_brdy.map_brdy.begin()->second.Nodet = { 1.0, 0.0 };	// tangent left - end  node
		(--bellhop_hydrology_out.top_brdy.map_brdy.end())->second.Nodet = { 1.0, 0.0 };// tangent right - end node

		for (auto iter_brdy = ++bellhop_hydrology_out.top_brdy.map_brdy.begin(); iter_brdy != --bellhop_hydrology_out.top_brdy.map_brdy.end(); iter_brdy++)
		{
			float sss = 0.5;	// странна вещь, но именнь так и написано
			iter_brdy->second.Nodet = (++iter_brdy)->second.t * sss + (--iter_brdy)->second.t * (1.0f - sss);
		}
		for (auto iter_brdy : bellhop_hydrology_out.top_brdy.map_brdy)
		{
			iter_brdy.second.Noden.r = iter_brdy.second.Nodet.z;
			iter_brdy.second.Noden.z = -iter_brdy.second.Nodet.r;
		}
		for (auto iter_brdy = bellhop_hydrology_out.top_brdy.map_brdy.begin(); iter_brdy != --bellhop_hydrology_out.top_brdy.map_brdy.end(); iter_brdy++)
		{
			// compute curvature in each segment
			auto phi1 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy++;
			auto phi2 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy--;
			iter_brdy->second.Kappa = (phi2 - phi1) / iter_brdy->second.Len;	// this is curvature = dphi / ds
		}
	}
	else
	{
		for (auto iter_brdy : bellhop_hydrology_out.top_brdy.map_brdy)
		{
			iter_brdy.second.Kappa = 0;
		}
	}

	// дно
	bellhop_hydrology_out.bot_brdy.type = this->bot_brdy.type;
	for (auto iter_r : this->bot_brdy.map_brdy)
	{
		auto _brdy = iter_r.second;
		_brdy.x.r = _r - iter_r.first;
		bellhop_hydrology_out.bot_brdy.map_brdy.insert({ _brdy.x.r, _brdy });
		
		if (iter_r.first > _r)
		{
			_brdy.x.r = 0;

			bellhop_hydrology_out.bot_brdy.map_brdy.insert({ _brdy.x.r, _brdy });
		}
	}

	for (std::map<float, bdryMod::brdy>::iterator iter_brdy = this->bot_brdy.map_brdy.begin(); iter_brdy != --this->bot_brdy.map_brdy.end(); iter_brdy++)
	{
		coord_2d x = (++iter_brdy)->second.x;
		coord_2d y = (--iter_brdy)->second.x;
		iter_brdy->second.t = x - y;
		iter_brdy->second.Len = iter_brdy->second.t.abs();
		iter_brdy->second.t = iter_brdy->second.t / iter_brdy->second.Len;

		iter_brdy->second.normal.r = -iter_brdy->second.t.z;
		iter_brdy->second.normal.z = +iter_brdy->second.t.r;
	}

	if (bellhop_hydrology_out.bot_brdy.type == bdryMod::types::curve)	// curvilinear option : compute tangent and normal at node by averaging normals on adjacent segments
	{
		bellhop_hydrology_out.bot_brdy.map_brdy.begin()->second.Nodet = { 1.0, 0.0 };	// tangent left - end  node
		(--bellhop_hydrology_out.bot_brdy.map_brdy.end())->second.Nodet = { 1.0, 0.0 };// tangent right - end node

		for (auto iter_brdy = ++bellhop_hydrology_out.bot_brdy.map_brdy.begin(); iter_brdy != --bellhop_hydrology_out.bot_brdy.map_brdy.end(); iter_brdy++)
		{
			float sss;
			//			sss = (--iter_brdy)->second.Len / (iter_brdy->second.Len + (++iter_brdy)->second.Len);
			sss = 0.5;	// странна вещь, но именнь так и написано
			iter_brdy->second.Nodet = (++iter_brdy)->second.t * sss + (--iter_brdy)->second.t * (1.0f - sss);
		}
		for (auto iter_brdy : bellhop_hydrology_out.bot_brdy.map_brdy)
		{
			iter_brdy.second.Noden.r = -iter_brdy.second.Nodet.z;
			iter_brdy.second.Noden.z = iter_brdy.second.Nodet.r;
		}
		for (auto iter_brdy = bellhop_hydrology_out.bot_brdy.map_brdy.begin(); iter_brdy != --bellhop_hydrology_out.bot_brdy.map_brdy.end(); iter_brdy++)
		{
			// compute curvature in each segment
			auto phi1 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy++;
			auto phi2 = atan2(iter_brdy->second.Nodet.z, iter_brdy->second.Nodet.r);
			iter_brdy--;
			iter_brdy->second.Kappa = (phi2 - phi1) / iter_brdy->second.Len;	// this is curvature = dphi / ds
		}
	}
	else
	{
		for (auto iter_brdy : bellhop_hydrology_out.bot_brdy.map_brdy)
		{
			iter_brdy.second.Kappa = 0;
		}
	}

	return bellhop_hydrology_out;
};
