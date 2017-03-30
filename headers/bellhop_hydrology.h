#ifndef HAC_HYDROLOGY_H
#define HAC_HYDROLOGY_H

#include "../headers/main.h"
#include "../headers/ssp.h"
#include "../headers/brdy.h"
#include "../headers/bands.h"

#include "../../GISDataBase/headers/det_hydro.h"

HAC_LIBRARY_NAMESPACE
{
	struct HAC_EXPORT bellhop_hydrology	//! структура для описания гидрологии в библиотеке НАС
	{
		unsigned int heavy_sea;	//!< волнение моря, баллы
		freq_band FB;	//!< частота, Гц
		float koef_attenuation;	//!< коэффициент затухания, дБ/км		
		bdryMod top_brdy;	//!< поверхность
		bdryMod bot_brdy;	//!< дно
		sspmod ssp;	//!< набор ВРСЗ
		std::map<float, std::map<float, float>> scat;	//!< 1 - расстояние м, 2 - глубина м, 3 - коэффициент объемного рассеяния на запрашиваемой глубине 1/м^3

		/*!
		\param[in] _hydrology задаваемая из базы гидрология
		\param[in] _ati_type тип интерполяции для поверхности моря
		\param[in] _bty_type тип интерполяции для дна моря
		*/
		error const set_param(const GDB::det::hydrology &_hydrology, const bdryMod::types &_ati_type, const bdryMod::types &_bty_type);	//!< задание гидрологии из базы

		bellhop_hydrology reverse(const float &_r);	// для помехи нужна гидрология в обратном порядке (_r - расстояние от предыдущего нуля до точки расчета помехи)
	};
};

#endif /*HAC_HYDROLOGY_H*/